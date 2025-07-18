# Config wrappers in LUMI

An interactive shell session inside the container can be started with the
following command:

```bash
$ ./config/csc/exec.lumi.container.cce bash
```

where `exec.lumi.container.cce` is a symbolic link pointing to the current
default container wrapper.

Once the shell session is started, ICON can be configured and built as follows:

- for CPU:
    ```bash
    Apptainer> ./config/csc/lumi.cpu.cce
    ...
    Apptainer> make -j8
    ```
- for GPU:
    ```bash
    Apptainer> ./config/csc/lumi.gpu.cce
    ...
    Apptainer> make -j8
    ```

where `lumi.cpu.cce` and `lumi.gpu.cce` are symbolic links pointing the current
default configure wrappers.

Alternatively, ICON can be configured, built and run as usual in the native
environment using the `lumi.*.native.*` configure wrappers.
