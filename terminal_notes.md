**Counting the number of genes across fasta**

Each gene starts with a specific ">" charecter. To count the number of 
instances of this charecter in the file, I use:

``` bash
grep -c ">" <file.name>
```


