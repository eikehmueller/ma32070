Start interactive session with

```
srun --partition spot-hbv3-120 --nodes 1  --account VB-MA5HPC-011  --qos spot-hbv3-120 --job-name "interactive" --cpus-per-task 120 --time 01:00:00 --pty bash
```

Load modules with
```
source load_modules.bash
```