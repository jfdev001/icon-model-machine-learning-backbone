- [ ] I have provided the URL to the commit used for the failed build.
- [ ] I have attached the log files.

    You can generate a tarball with the log files by running the following commands from the root build directory of ICON (will generate a `build-report.tar.gz`):
    ```bash
    make V=1 2>&1 | tee make.log
    tar --transform 's:^:build-report/:' -czf build-report.tar.gz $(find . -name 'config.log') make.log
    ```
