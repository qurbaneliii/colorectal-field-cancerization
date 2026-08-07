# Docker verification

Status: **BLOCKED**

The Docker 29.5.3 client and Docker Desktop are installed. A fresh `docker info`
check on 2026-08-08 returned exit code 1 because the Linux-engine named pipe
`//./pipe/dockerDesktopLinuxEngine` does not exist. On 2026-08-01 all
three required commands were attempted and returned exit code 1:

- `docker info`
- `docker build -t colorectal-field-cancerization:publication .`
- `docker run --rm colorectal-field-cancerization:publication make smoke`

Each failed before execution because the Linux-engine named pipe
`//./pipe/dockerDesktopLinuxEngine` does not exist. The earlier host diagnostic
`wsl.exe --status` returned `Wsl/0x80070422`, indicating the required Windows
service is disabled or has no enabled device associated with it.

No image-build success is claimed. An administrator must enable the host WSL
service/feature and start Docker Desktop, then run:

```powershell
docker info
docker build -t colorectal-field-cancerization:publication .
docker run --rm colorectal-field-cancerization:publication make smoke
```

Replace this report with the server version, image digest, build exit code, and
container verification output only after those commands succeed.
