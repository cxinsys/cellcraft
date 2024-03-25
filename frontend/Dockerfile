# 1단계: 빌드 환경 설정
# Node.js 이미지를 기반으로 설정합니다.
FROM node:lts-alpine as build-stage

# 작업 디렉토리를 설정합니다.
WORKDIR /app

# 애플리케이션 의존성 파일들을 컨테이너로 복사합니다.
COPY package*.json ./

# 의존성을 설치합니다.
RUN npm install

# 애플리케이션 소스 코드를 컨테이너로 복사합니다.
COPY . .

# 애플리케이션을 빌드하여 정적 파일을 생성합니다.
RUN npm run build
