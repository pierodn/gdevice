BUILD.md

Purpose
-------
Short instructions to build this project on Windows (Visual Studio / MSBuild).

Prerequisites
-------------
- Windows 10/11
- Visual Studio (2008+) with C++ workload or MSBuild + Windows SDK

Quick build (Developer Command Prompt or PowerShell)
--------------------------------------------------
1. Build with MSBuild (Release, Win32):
   msbuild /m /p:Configuration=Release /p:Platform=Win32 .\gdevice.2008.sln

Notes
-----
- The solution targets Win32; ensure Platform=Win32 when building.
- Tests/walker consumes assets under tests/walker/assets and expects to be run from repository root or with working directory set to tests\walker.
- A helper tool (bin2c) is included as a project; build it first if you regenerate embedded data.
- Repository currently contains committed binaries in /binaries. Consider removing them from VCS (use .gitignore and Git LFS or Releases) to reduce repo size.

CI
--
The GitHub Actions workflow (.github/workflows/msbuild.yml) currently sets SOLUTION_FILE_PATH to ./source/gdevice.2008.sln; change it to ./gdevice.2008.sln to match this repo.

Next recommended steps
----------------------
1. Add a .gitignore and remove committed binaries (or migrate to Git LFS).
2. Update CI workflow path and add artifact upload step.
3. Add a short README entry referencing this BUILD.md.