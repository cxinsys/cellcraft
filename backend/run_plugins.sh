#!/bin/bash
# run_plugins.sh

# plugins 폴더 내의 모든 .sh 파일에 대해 실행
for script in /app/plugins/*.sh; do
    chmod +x "$script"  # 실행 권한 부여
    "$script"           # 스크립트 실행
done
