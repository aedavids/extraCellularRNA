# Need to move our salom index. 
The web based file upload does not work for large files or if you have a large number of files

ref: [moving data to an and from google bucket](https://support.terra.bio/hc/en-us/articles/360024056512-Moving-data-to-from-a-workspace-or-external-Google-bucket-#h_01EN30QXSJ6HEYN4GYKDFR9Q8D)

I found the project name on the tera workspace dashboard page
```
(extraCellularRNA) $ gsutil ls -p test-aedavids-proj
gs://cromwell-auth-test-aedavids-proj/
gs://fc-6c27c400-9457-44a7-bcac-83c5c8228857/
gs://fc-eaaf45fa-ba04-4361-9ae5-4f70dc8bbab8/
gs://fc-f77cf2c4-3bef-47ac-8d19-3c43cc61574b/
gs://fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059/
gs://storage-logs-test-aedavids-proj/
```

Google bucket. I found this on the terra workspace dashboard page
fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059

```
(extraCellularRNA) $ gsutil ls gs://fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059
gs://fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059/git.commit.msg.txt
gs://fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059/6a6c9b92-3026-47d3-8944-60f0842c566e/
gs://fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059/af96753d-ff8d-4445-915e-55f97543b137/
gs://fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059/d44c8949-6909-4fe9-85e1-e88d2143ca51/
```

Uploading large files using parallel composite uploads
TODO: https://cloud.google.com/storage/docs/uploads-downloads#parallel-composite-uploads
```
(extraCellularRNA) $ gsutil cp sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.tar.gz gs://fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059
Copying file://sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.tar.gz [Content-Type=application/x-tar]...
If you experience problems with multiprocessing on MacOS, they might be related to https://bugs.python.org/issue33725. You can disable multiprocessing by editing your .boto config or by adding the following flag to your command: `-o "GSUtil:parallel_process_count=1"`. Note that multithreading is still available even if you disable multiprocessing.

If you experience problems with multiprocessing on MacOS, they might be related to https://bugs.python.org/issue33725. You can disable multiprocessing by editing your .boto config or by adding the following flag to your command: `-o "GSUtil:parallel_process_count=1"`. Note that multithreading is still available even if you disable multiprocessing.

/ [1 files][ 16.0 GiB/ 16.0 GiB]  375.3 KiB/s                                   
Operation completed over 1 objects/16.0 GiB.                                     
(extraCellularRNA) $ gsutil cp sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.tar.gz gs://fc-secure-519db2bc-049f-43a0-ab75-a2eb9c2cb059
Copying file://sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.tar.gz [Content-Type=application/x-tar]...
If you experience problems with multiprocessing on MacOS, they might be related to https://bugs.python.org/issue33725. You can disable multiprocessing by editing your .boto config or by adding the following flag to your command: `-o "GSUtil:parallel_process_count=1"`. Note that multithreading is still available even if you disable multiprocessing.

If you experience problems with multiprocessing on MacOS, they might be related to https://bugs.python.org/issue33725. You can disable multiprocessing by editing your .boto config or by adding the following flag to your command: `-o "GSUtil:parallel_process_count=1"`. Note that multithreading is still available even if you disable multiprocessing.

- [1 files][ 16.0 GiB/ 16.0 GiB]  369.1 KiB/s                                   
Operation completed over 1 objects/16.0 GiB.                                     
(extraCellularRNA) $ 

```
