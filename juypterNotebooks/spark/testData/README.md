
the model files where crafted by hand. The trick to over ridding the tab setting in .emacs file is

in the emacs scratch buffer. copy the following line

```
(setq-default indent-tabs-mode 'only)
```

Move the cursor to the end of the line

type ctrl-j

You can also use tr to convert white space into tab

```
$ cat sample.tsv | tr  --squeeze-repeats ' ' '\t' > f
```

# sparkDESeqTest2 contains a more complete set of test files we can use to debug our spark code eqivalent of DESEq2, tximport, ...
