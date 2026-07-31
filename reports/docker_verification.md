# Docker verification

Status: **BLOCKED**

The Docker 29.5.3 client and Docker Desktop are installed. `docker info` could
not reach the Linux engine. After requesting a Docker Desktop start,
`docker desktop status` remained `starting`, and `wsl.exe --status` returned
`Wsl/0x80070422`: the required Windows service is disabled or has no enabled
device associated with it.

No image-build success is claimed. An administrator must enable the host WSL
service/feature and start Docker Desktop, then run:

```powershell
docker info
docker build -t colorectal-field-cancerization:publication .
docker run --rm colorectal-field-cancerization:publication
```

Replace this report with the server version, image digest, build exit code, and
container verification output only after those commands succeed.
