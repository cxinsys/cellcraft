FROM tiangolo/uvicorn-gunicorn:python3.10

# 필요한 패키지 설치
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    git \
    build-essential \
    autoconf \
    automake \
    libtool \
    curl \
    make \
    g++ \
    unzip \
    default-jdk \
    python3-pip \
    python3-dev \
    jq \
    && rm -rf /var/lib/apt/lists/*

# Miniconda 설치
RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm Miniconda3-latest-Linux-x86_64.sh

# Conda 환경 변수 설정
ENV PATH /opt/conda/bin:$PATH

# Conda 환경 설정
COPY environment.yml ./
RUN conda env create -f environment.yml

# 작업 디렉토리 설정
WORKDIR /app

# 애플리케이션 소스 코드 복사
COPY . .

# Conda 환경 활성화 및 FastAPI 애플리케이션 실행
CMD ["/opt/conda/envs/snakemake/bin/python", "-m", "uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
