**Counting the number of genes across fasta**

Each gene starts with a specific ">" character. To count the number of 
instances of this charecter in the file, I use:

``` bash
grep -c ">" <file.name>
```

**Terminal and session control**

`screen`:
``` bash 
screen -S <name> # create the screen session with the <name>
screen -r <name>/<session number> # reattach to the sesssion with <name>
screen -ls # list available screen sessions 
```
Ctrl+a+d -> detach from the screen


Check all the processes in the session:
``` bash 
ps # list all processes currently running under your user in the session
ps aux | grep <script name> # check specific process or script


jobs # list running jobs
htop # check CPU usage
```



