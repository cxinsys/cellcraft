# cellcraft
> bio web Backend

# install dependencies
pip install -r requirements.txt

# virtual environment settings
python -m venv .venv

# activate venv
source venv/bin/activate

# uvicorn server with reload at 127.0.0.1:8000
uvicorn app.main:app --reload