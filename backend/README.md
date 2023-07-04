# cellcraft
> bio web Backend

## Build Setup


1. install anaconda virtual environment
``` bash
conda env create -f snakemake.yaml
conda activate snakemake
```
2. install dependencies
``` bash
pip install -r requirements.txt
```
3. DB (PostgreSQL connection)

[Install PostgreSQL](https://www.postgresql.org/download/)
``` bash
psql postgres
create database cellcraft;
```
4. uvicorn server at 127.0.0.1:8000
``` bash
uvicorn app.main:app
```
5. rabbitMQ server install & start
```bash
brew install rabbitmq
brew services start rabbitmq
```
6. celery worker run
```bash
celery -A app.main.celery worker --loglevel=info -E -Q workflow_task
```
