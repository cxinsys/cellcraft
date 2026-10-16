# cellcraft
> bio web Backend

## Build Setup


1. install anaconda virtual environment
``` bash
conda env create -f environment.yml
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
4. environment variables (`backend/.env`, read by `app/core/config.py`)

Defaults point at docker-compose service names (`rabbitmq`, `redis`), so override them for a local run.
```bash
POSTGRES_USER=<user>
POSTGRES_PASSWORD=<password>
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=cellcraft
CELERY_BROKER_URL=amqp://guest:guest@localhost:5672//
REDIS_URL=redis://localhost:6379/0
```
5. uvicorn server at 127.0.0.1:8000
``` bash
uvicorn app.main:app
```
6. rabbitMQ server install & start
```bash
brew install rabbitmq
brew services start rabbitmq
```
7. celery worker run
```bash
celery -A app.worker.app worker --loglevel=info -E -Q workflow_task,plugin_task
```
