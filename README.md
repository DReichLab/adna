aDNA is a set of tools to process ancient DNA data, in particular data produced
with the [ReichLab protocol][udg]. As of now, aDNA only consists of one tool,
adna-trim. It efficiently trims sequencing adapters, checks inline barcodes and
merges overlapping ends, all in one go. The typical command line to invoke
adna-trim is:
```sh
seqtk mergepe R1.fq.gz R2.fq.gz | adna-trim -p out-pe -b barcode.txt - | gzip -1 > out-se.fq.gz
```
where [seqtk][seqtk] generates an interleaved FASTQ and `barcode.txt` gives the
read1 and read2 barcodes (see [example.bc][example-bc] for an example). Read
pairs that can't be unambiguously merged will be written to `out-pe.R1.fq.gz`
and `out-pe.R2.fq.gz`; merged reads will be outputted to the standard output.
Inline barcodes are appended to read names.

[udg]: http://rstb.royalsocietypublishing.org/content/370/1660/20130624
[seqtk]: https://github.com/lh3/seqtk
[example-bc]: https://github.com/DReichLab/adna/blob/master/example.bc
