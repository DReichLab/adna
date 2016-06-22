#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include <stdio.h>
#include "kvec.h"
#include "kstring.h"

/********************
 * Global variables *
 ********************/

const char *lt_adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"; // Illumina 3'-end adapter
const char *lt_adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";

#define AD_MAX_BC_LEN 16
#define AD_N_BARCODE 4

const char *ad_barcode_for[2][AD_N_BARCODE] = { { "CTAGACA", "GACTCGC", "TCGAGTG", "AGTCTAT" }, { "CCTGCAG", "GGATGCT", "TTCATGA", "AAGCATC" } };
const char *ad_barcode_rev[2][AD_N_BARCODE] = { { "TGTCTAG", "GCGAGTC", "CACTCGA", "ATAGACT" }, { "CTGCAGG", "AGCATCC", "TCATGAA", "GATGCTT" } };

#define AD_MAX_TYPE 30

enum ad_type_e {
	AD_UNKNOWN = 0,
	AD_NO_BARCODE = 1,
	AD_MAL_ADAP = 2,
	AD_NO_ADAP = 3,
	AD_SHORT_SE = 11,
	AD_COMPLETE_MERGE = 21,
	AD_PARTIAL_MERGE = 22,
	AD_AMBI_MERGE = 23,
	AD_NO_MERGE = 24
};

typedef struct {
	int n_threads;
	int chunk_size;
	int min_seq_len;
	int max_qual;
	int max_ovlp_pen, min_ovlp_len;
	int max_adap_pen, min_adap_len;
	int bc_len;
	int tab_out;
} lt_opt_t;

static void lt_opt_init(lt_opt_t *opt)
{
	memset(opt, 0, sizeof(lt_opt_t));
	opt->n_threads = 2;
	opt->chunk_size = 10000000;
	opt->max_qual = 50;
	opt->min_seq_len = 30;
	opt->max_ovlp_pen = 2;
	opt->min_ovlp_len = 10;
	opt->max_adap_pen = 1;
	opt->min_adap_len = 3;
	opt->bc_len = 7;
}

/******************
 * K-mer matching *
 ******************/

#include "khash.h"
KHASH_SET_INIT_INT64(s64)
typedef khash_t(s64) lt_seqcloud1_t;

typedef struct {
	int l;
	uint64_t s;
	lt_seqcloud1_t *mm;
} lt_seqcloud_t;

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static lt_seqcloud_t *lt_sc_init(void)
{
	lt_seqcloud_t *sc;
	sc = (lt_seqcloud_t*)calloc(1, sizeof(lt_seqcloud_t));
	sc->mm = kh_init(s64);
	return sc;
}

void lt_sc_destroy(lt_seqcloud_t *sc)
{
	kh_destroy(s64, sc->mm);
	free(sc);
}

void lt_sc_add_core(lt_seqcloud_t *sc, uint64_t s)
{
	int i, absent;
	sc->s = s = s & ((1ULL<<sc->l*2) - 1);
	for (i = 0; i < sc->l; ++i) {
		int i2 = i * 2, a, c = s>>i2&3;
		for (a = 1; a < 4; ++a) {
			uint64_t x = (s & ~(3ULL << i2)) | (uint64_t)((a+c)&3) << i2;
			kh_put(s64, sc->mm, x, &absent);
		}
	}
}

lt_seqcloud_t *lt_sc_gen(const char *s)
{
	lt_seqcloud_t *sc;
	uint64_t x = 0;
	int i;
	sc = lt_sc_init();
	sc->l = strlen(s);
	for (i = 0; s[i] && i < sc->l; ++i) {
		int c = seq_nt4_table[(uint8_t)s[i]];
		if (c > 3) {
			lt_sc_destroy(sc);
			return 0;
		}
		x = x << 2 | c;
	}
	lt_sc_add_core(sc, x);
	return sc;
}

typedef struct {
	uint32_t pos:30, type:2;
} lt_sc_hit_t;

int lt_sc_test(const lt_seqcloud_t *sc, const char *seq, int max_hits, lt_sc_hit_t *hits)
{
	int i, l, n = 0;
	uint64_t x = 0, mask = (1ULL << sc->l*2) - 1;
	for (i = l = 0; seq[i]; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x = (x << 2 | c) & mask;
			if (++l >= sc->l) {
				if (x == sc->s || kh_get(s64, sc->mm, x) != kh_end(sc->mm)) {
					hits[n].pos = i - (sc->l - 1);
					hits[n++].type = x == sc->s? 0 : 1;
					if (n == max_hits) return n;
				}
			}
		} else l = 0, x = 0;
	}
	return n;
}

/**********************
 * Reverse complement *
 **********************/

char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

void lt_seq_rev(int l, const char *f, char *r)
{
	int i;
	for (i = 0; i < l; ++i)
		r[l - i - 1] = f[i];
	r[l] = 0;
}

void lt_seq_revcomp(int l, const char *f, char *r)
{
	int i;
	for (i = 0; i < l; ++i)
		r[l - i - 1] = (uint8_t)f[i] >= 128? 'N' : comp_tab[(uint8_t)f[i]];
	r[l] = 0;
}

/**********************
 * Ungapped extension *
 **********************/

#define LT_QUAL_THRES 53 // =33+20
#define LT_HIGH_PEN 3
#define LT_LOW_PEN  1

int lt_ue_for1(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int min_len, int max_pen)
{
	int i, pen = 0;
	for (i = 0; i < l1 && i < l2; ++i) {
		if (s1[i] != s2[i]) {
			pen += q1[i] >= LT_QUAL_THRES && (q2 == 0 || q2[i] >= LT_QUAL_THRES)? LT_HIGH_PEN : LT_LOW_PEN;
			if (i <= min_len && pen > max_pen) break;
			if (i > min_len && pen * min_len > i * max_pen) break; // in effect: pen > max_pen * ((double)i / min_len)
		}
	}
	return i;
}

int lt_ue_rev1(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int min_len, int max_pen)
{
	int i, pen = 0;
	for (i = 0; i < l1 && i < l2; ++i) {
		if (s1[l1-1-i] != s2[l2-1-i]) {
			pen += q1[l1-1-i] >= LT_QUAL_THRES && (q2 == 0 || q2[l2-1-i] >= LT_QUAL_THRES)? LT_HIGH_PEN : LT_LOW_PEN;
			if (i <= min_len && pen > max_pen) break;
			if (i > min_len && pen * min_len > i * max_pen) break;
		}
	}
	return i;
}

int lt_ue_for(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int max_pen, int min_len, int max_pos, uint64_t *pos)
{
	int i, n = 0;
	for (i = min_len; i <= l1; ++i) {
		int l;
		l = lt_ue_for1(i, s1 + l1 - i, q1 + l1 - i, l2, s2, q2, min_len, max_pen);
		if (l >= min_len && (l == i || l == l2)) {
			pos[n++] = (uint64_t)(l1 - i) << 32 | l;
			if (n == max_pos) return n;
		}
	}
	return n;
}

int lt_ue_rev(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int max_pen, int min_len, int max_pos, uint64_t *pos)
{
	int i, n = 0;
	for (i = min_len; i <= l1; ++i) {
		int l;
		l = lt_ue_rev1(i, s1, q1, l2, s2, q2, min_len, max_pen);
		if (l >= min_len && (l == i || l == l2)) {
			pos[n++] = (uint64_t)(l1 - i) << 32 | l;
			if (n == max_pos) return n;
		}
	}
	return n;
}

int lt_ue_contained(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int max_pen, int max_pos, uint64_t *pos)
{
	int i, n = 0;
	for (i = 1; i < l2 - l1; ++i) {
		int l;
		l = lt_ue_for1(l1, s1, q1, l2 - i, s2 + i, q2 + i, l1, max_pen);
		if (l == l1) {
			pos[n++] = (uint64_t)i << 32 | l;
			if (n == max_pos) return n;
		}
	}
	return n;
}

/**********************
 * Batch FASTQ reader *
 **********************/

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	uint32_t l_seq:31, dbl_bind:1;
	enum ad_type_e type;
	char *name, *seq, *qual, *bc;
} bseq1_t;

bseq1_t *bseq_read(kseq_t *ks, int chunk_size, int *n_)
{
	int size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->seq = strdup(ks->seq.s);
		s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		s->bc = 0;
		s->l_seq = ks->seq.l;
		s->dbl_bind = 0;
		s->type = AD_UNKNOWN;
		size += seqs[n++].l_seq;
		if (size >= chunk_size && (n&1) == 0) break;
	}
	*n_ = n;
	return seqs;
}

/*********************************
 * Core trimming/merging routine *
 *********************************/

typedef struct {
	lt_opt_t opt;
	lt_seqcloud_t *bc_for[2][AD_N_BARCODE], *bc_rev[2][AD_N_BARCODE];
	kseq_t *ks;
	gzFile fp_pe[2];
	uint64_t types[AD_MAX_TYPE+1];
} lt_global_t;

void lt_global_init(lt_global_t *g)
{
	int k, i;
	memset(g, 0, sizeof(lt_global_t));
	lt_opt_init(&g->opt);
	for (k = 0; k < 2; ++k) {
		for (i = 0; i < AD_N_BARCODE; ++i) {
			g->bc_for[k][i] = lt_sc_gen(ad_barcode_for[k][i]);
			g->bc_rev[k][i] = lt_sc_gen(ad_barcode_rev[k][i]);
		}
	}
}

#define MAX_BINDING_HITS 3

static inline void trim_bseq_5(bseq1_t *s, int l)
{
	memmove(s->seq, s->seq + l, s->l_seq - l);
	memmove(s->qual, s->qual + l, s->l_seq - l);
	s->l_seq -= l;
	s->seq[s->l_seq] = s->qual[s->l_seq] = 0;
}

static inline int merge_base(int max_qual, char fc, char fq, char rc, char rq)
{
	int y;
	if (fc == rc) {
		int q = fq > rq? fq - 33 : rq - 33;
		y = toupper(fc) | (33 + (q < max_qual? q : max_qual)) << 8;
	} else {
		if (fq > rq) y = toupper(fc) | (33 + (fq - rq)) << 8;
		else y = toupper(rc) | (33 + (rq - fq)) << 8;
	}
	return y;
}

static inline void trim_adap(bseq1_t *s, const char *adap, int is_5, int min_len, int max_pen, int allow_contained)
{
	int n_hits, l_adap;
	uint64_t hits[4];
	l_adap = strlen(adap);
	if (is_5) n_hits = lt_ue_rev(s->l_seq, s->seq, s->qual, l_adap, adap, 0, max_pen, min_len, 4, hits);
	else n_hits = lt_ue_for(s->l_seq, s->seq, s->qual, l_adap, adap, 0, max_pen, min_len, 4, hits);
	if (n_hits > 0 && (allow_contained || (hits[0]>>32) + (uint32_t)hits[0] == s->l_seq || (hits[n_hits-1]>>32) + (uint32_t)hits[n_hits-1] == s->l_seq)) {
		int len = s->l_seq - (hits[n_hits-1]>>32); // trim the longest hit
		if (is_5) {
			if (len > min_len) trim_bseq_5(s, len); // trim
			else memset(s->qual, 33+1, len); // reduce baseQ
		} else {
			if (len > min_len) s->l_seq -= len, s->seq[s->l_seq] = s->qual[s->l_seq] = 0; // trim
			else memset(s->qual + (s->l_seq - len), 33+1, len); // reduce baseQ
		}
	}
}

void lt_process(const lt_global_t *g, bseq1_t s[2])
{
	int i, k, mlen, olen[2], bc_id[2];
	char *rseq, *rqual, *xseq, *xqual, *bc;

	mlen = s[0].l_seq > s[1].l_seq? s[0].l_seq : s[1].l_seq;
	rseq = (char*)alloca(mlen + 1);
	rqual = (char*)alloca(mlen + 1);
	xseq = (char*)alloca(s[0].l_seq + s[1].l_seq + 1);
	xqual = (char*)alloca(s[0].l_seq + s[1].l_seq + 1);
	bc = (char*)alloca(mlen + 1);

	// trim trailing N
	for (k = 0; k < 2; ++k) {
		bseq1_t *sk = &s[k];
		for (i = sk->l_seq - 1; i >= 0; --i) // trim trailing "N"
			if (sk->seq[i] != 'N') break;
		sk->l_seq = i + 1;
		sk->seq[sk->l_seq] = sk->qual[sk->l_seq] = 0;
		olen[k] = sk->l_seq;
	}
	// test barcode
	for (k = 0; k < 2; ++k) {
		char bc[AD_MAX_BC_LEN + 1];
		lt_sc_hit_t hits[3];
		strncpy(bc, s[k].seq, g->opt.bc_len);
		bc[g->opt.bc_len] = 0;
		for (i = 0; i < AD_N_BARCODE; ++i)
			if (lt_sc_test(g->bc_for[k][i], bc, 3, hits))
				break;
		if (i == AD_N_BARCODE) break;
		bc_id[k] = i;
	}
	if (k < 2) {
		s[0].type = s[1].type = AD_NO_BARCODE;
		return;
	}
	// write barcode
	s->bc = (char*)calloc(g->opt.bc_len * 2 + 1, 1);
	strncpy(s->bc, ad_barcode_for[0][bc_id[0]], g->opt.bc_len);
	strncpy(s->bc + g->opt.bc_len, ad_barcode_for[1][bc_id[1]], g->opt.bc_len);
	trim_bseq_5(&s[0], g->opt.bc_len);
	trim_bseq_5(&s[1], g->opt.bc_len);
	// trim Illumina PE adapters
	olen[0] = s[0].l_seq, olen[1] = s[1].l_seq;
	trim_adap(&s[0], lt_adapter1, 0, g->opt.min_adap_len, g->opt.max_adap_pen, 1);
	trim_adap(&s[1], lt_adapter2, 0, g->opt.min_adap_len, g->opt.max_adap_pen, 1);
	if (s[0].l_seq == olen[0] && s[1].l_seq == olen[1]) {
		s[0].type = s[1].type = AD_NO_ADAP;
	} else if (s[0].l_seq == s[1].l_seq && s[0].l_seq < olen[0] && s[1].l_seq < olen[1] && s->l_seq > g->opt.bc_len) {
		for (k = 0; k < 2; ++k) {
			s[k].type = AD_COMPLETE_MERGE;
			s[k].l_seq -= g->opt.bc_len;
			s[k].seq[s[k].l_seq] = s[k].qual[s[k].l_seq] = 0;
		}
		for (i = 0; i < s->l_seq; ++i) {
			int j = s->l_seq - 1 - i;
			int y, r = (uint8_t)s[1].seq[j] >= 128? 'N' : comp_tab[(uint8_t)s[1].seq[j]];
			y = merge_base(g->opt.max_qual, s[0].seq[i], s[0].qual[i], r, s[1].qual[j]);
			s->seq[i] = y & 0xff;
			s->qual[i] = y >> 8;
		}
		s[1].l_seq = 0;
		if (s->l_seq < g->opt.min_seq_len) s[0].type = s[1].type = AD_SHORT_SE;
	} else if (s[0].l_seq != s[1].l_seq || s[0].l_seq == olen[0] || s[1].l_seq == olen[1]) {
		s[0].type = s[1].type = AD_MAL_ADAP;
	}
	// find end overlaps
	if (s->type == AD_NO_ADAP) {
		int n_fh;
		uint64_t fh[2];
		// reverse the other read
		lt_seq_revcomp(s[1].l_seq, s[1].seq, rseq);
		lt_seq_rev(s[1].l_seq, s[1].qual, rqual);
		// find overlaps
		n_fh = lt_ue_for(s[0].l_seq, &s[0].seq[0], &s[0].qual[0], s[1].l_seq, rseq, rqual, g->opt.max_ovlp_pen, g->opt.min_ovlp_len, 2, fh);
		if (n_fh > 1) {
			s[0].type = s[1].type = AD_AMBI_MERGE;
		} else if (n_fh == 0) {
			s[0].type = s[1].type = AD_NO_MERGE;
		} else {
			int x = 0;
			s[0].type = s[1].type = AD_PARTIAL_MERGE;
			if (n_fh == 1) {
				int l = (uint32_t)fh[0], st = fh[0]>>32;
				for (i = 0; i < st; ++i)
					xseq[x] = s[0].seq[i], xqual[x++] = s[0].qual[i];
				for (i = 0; i < l; ++i) {
					int j = st + i, y;
					y = merge_base(g->opt.max_qual, s[0].seq[j], s[0].qual[j], rseq[i], rqual[i]);
					xseq[x] = (uint8_t)y, xqual[x++] = y>>8;
				}
				if (l < s[1].l_seq) {
					for (i = l; i < s[1].l_seq; ++i)
						xseq[x] = rseq[i], xqual[x++] = rqual[i];
				} else {
					for (i = l; i < s[0].l_seq; ++i)
						xseq[x] = s[0].seq[i], xqual[x++] = s[0].qual[i];
				}
			}
			xseq[x] = xqual[x] = 0;
			if (x < g->opt.min_seq_len) s[0].type = s[1].type = AD_SHORT_SE;
			free(s[0].seq); free(s[0].qual);
			s[0].seq = strdup(xseq);
			s[0].qual = strdup(xqual);
			s[0].l_seq = x;
			s[1].l_seq = 0;
		}
	}
	if (s->type == AD_NO_MERGE || s->type == AD_AMBI_MERGE)
		s[1].bc = strdup(s->bc);
}

/**********************
 * Callback functions *
 **********************/

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	lt_global_t *g;
} data_for_t;

static void worker_for(void *_data, long i, int tid)
{
	data_for_t *data = (data_for_t*)_data;
	lt_process(data->g, &data->seqs[i<<1]);
}

static void *worker_pipeline(void *shared, int step, void *_data)
{
	int i;
	lt_global_t *g = (lt_global_t*)shared;
	if (step == 0) {
		data_for_t *ret;
		ret = calloc(1, sizeof(data_for_t));
		ret->seqs = bseq_read(g->ks, g->opt.chunk_size, &ret->n_seqs);
		assert((ret->n_seqs&1) == 0);
		ret->g = g;
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		data_for_t *data = (data_for_t*)_data;
		kt_for(g->opt.n_threads, worker_for, data, data->n_seqs>>1);
		return data;
	} else if (step == 2) {
		data_for_t *data = (data_for_t*)_data;
		if (g->opt.tab_out) { // tabular output
			for (i = 0; i < data->n_seqs; i += 2) {
				bseq1_t *s = &data->seqs[i];
				printf("%s\t%d\n", s->name, s->type);
			}
		} else { // FASTQ output (FASTA not supported yet)
			kstring_t pe[2] = {{0,0,0},{0,0,0}};
			for (i = 0; i < data->n_seqs; i += 2)
				++g->types[data->seqs[i].type];
			for (i = 0; i < data->n_seqs; ++i) {
				bseq1_t *s = &data->seqs[i];
				if (s->l_seq <= 0) continue;
				if (g->fp_pe[0] && g->fp_pe[1]) {
					if (s->type == AD_PARTIAL_MERGE || s->type == AD_COMPLETE_MERGE) {
						printf("%c%s:%d:%s\n%s\n", s->qual? '@' : '>', s->name, s->type, s->bc, s->seq);
						if (s->qual) { puts("+"); puts(s->qual); }
					} else if (s->type == AD_NO_MERGE || s->type == AD_AMBI_MERGE) {
						kstring_t *str = &pe[i&1];
						ksprintf(str, "%c%s:%d:%s\n%s\n", s->qual? '@' : '>', s->name, s->type, s->bc, s->seq);
						if (s->qual) { kputs("+\n", str); kputs(s->qual, str); kputc('\n', str); }
					}
				} else {
					if (s->type == AD_PARTIAL_MERGE || s->type == AD_COMPLETE_MERGE || s->type == AD_NO_MERGE || s->type == AD_AMBI_MERGE) {
						printf("%c%s:%d:%s", s->qual? '@' : '>', s->name, s->type, s->bc);
						if (s->type != AD_PARTIAL_MERGE && s->type != AD_COMPLETE_MERGE) {
							putchar('/'); putchar("12"[i&1]);
						}
						putchar('\n');
						puts(s->seq);
						if (s->qual) { puts("+"); puts(s->qual); }
					}
				}
			}
			for (i = 0; i < 2; ++i) {
				gzwrite(g->fp_pe[i], pe[i].s, pe[i].l);
				free(pe[i].s);
			}
		}
		for (i = 0; i < data->n_seqs; ++i) { // deallocate
			bseq1_t *s = &data->seqs[i];
			free(s->bc); free(s->seq); free(s->qual); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

#include <unistd.h>

int main(int argc, char *argv[])
{
	int c;
	lt_global_t g;
	gzFile fp;
	char *pe_prefix = 0;

	lt_global_init(&g);
	while ((c = getopt(argc, argv, "Tt:b:l:o:p:")) >= 0) {
		if (c == 't') g.opt.n_threads = atoi(optarg);
		else if (c == 'T') g.opt.tab_out = 1;
		else if (c == 'b') g.opt.bc_len = atoi(optarg);
		else if (c == 'l') g.opt.min_seq_len = atoi(optarg);
		else if (c == 'o') g.opt.min_ovlp_len = atoi(optarg);
		else if (c == 'p') pe_prefix = optarg;
	}
	if (argc - optind < 1) {
		fprintf(stderr, "Usage: seqtk mergepe <read1.fq> <read2.fq> | adna-trim [options] -\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT     number of threads [%d]\n", g.opt.n_threads);
//		fprintf(stderr, "  -b INT     barcode length [%d]\n", g.opt.bc_len);
		fprintf(stderr, "  -l INT     min read/fragment length to output [%d]\n", g.opt.min_seq_len);
		fprintf(stderr, "  -o INT     min overlap length [%d]\n", g.opt.min_ovlp_len);
		fprintf(stderr, "  -p STR     output PE reads to STR.R[12].fq.gz [stdout]\n");
		fprintf(stderr, "  -T         tabular output for debugging\n");
		return 1;
	}

	// TODO: check the barcode length!!

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	g.ks = kseq_init(fp);
	if (pe_prefix) {
		kstring_t str = {0,0,0};
		for (c = 0; c < 2; ++c) {
			str.l = 0;
			ksprintf(&str, "%s.R%d.fq.gz", pe_prefix, c+1);
			g.fp_pe[c] = gzopen(str.s, "w1");
		}
	}

	kt_pipeline(2, worker_pipeline, &g, 3);
	for (c = 0; c <= AD_MAX_TYPE; ++c)
		if (g.types[c] > 0)
			fprintf(stderr, "[M::%s] %2d  %ld\n", __func__, c, (long)g.types[c]);

	if (g.fp_pe[0]) gzclose(g.fp_pe[0]);
	if (g.fp_pe[1]) gzclose(g.fp_pe[1]);
	kseq_destroy(g.ks);
	gzclose(fp);
	return 0;
}
