"""Operator CLI for on-premise installs (#138).

Run inside the backend container, from the backend root:

    docker compose -f <compose-file> exec backend \
        /opt/micromamba/envs/snakemake/bin/python -m app.cli reset-password --email <email>

Commands:
    create-admin    create an admin account (random password unless --password is given)
    reset-password  issue a temporary password; the user must change it at next login
"""
import argparse
import sys

from app.core.bootstrap import DEFAULT_ADMIN_USERNAME, generate_password, validate_new_password
from app.core.exceptions import ValidationFailedError
from app.core.security import get_password_hash
from app.db.base_registry import Base  # noqa: F401  (registers every model for relationships)
from app.db.session import get_db_session
from app.user import models


def create_admin(db, *, email: str, username: str, password: str = None) -> str:
    if db.query(models.User).filter_by(email=email).first():
        raise ValidationFailedError(f"Email already registered: {email}")
    generated = password is None
    if generated:
        password = generate_password()
    else:
        validate_new_password(password)
    db.add(models.User(
        username=username, email=email, hashed_password=get_password_hash(password),
        is_active=True, is_superuser=True, must_change_password=generated,
    ))
    db.commit()
    return password


def reset_password(db, *, email: str) -> str:
    user = db.query(models.User).filter_by(email=email).first()
    if not user:
        raise ValidationFailedError(f"User not found: {email}")
    password = generate_password()
    user.hashed_password = get_password_hash(password)
    user.must_change_password = True
    user.is_active = True
    db.commit()
    return password


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m app.cli", description="CellCraft operator commands")
    commands = parser.add_subparsers(dest="command", required=True)

    create = commands.add_parser("create-admin", help="create an admin account")
    create.add_argument("--email", required=True)
    create.add_argument("--username", default=DEFAULT_ADMIN_USERNAME)
    create.add_argument("--password", help="omit to generate a temporary password")

    reset = commands.add_parser("reset-password", help="issue a temporary password")
    reset.add_argument("--email", required=True)

    args = parser.parse_args(argv)
    try:
        with get_db_session() as db:
            if args.command == "create-admin":
                password = create_admin(db, email=args.email, username=args.username, password=args.password)
                print(f"Admin account created: {args.email}")
                if args.password is None:
                    print(f"Temporary password: {password}  (must be changed at first login)")
            else:
                password = reset_password(db, email=args.email)
                print(f"Temporary password for {args.email}: {password}  (must be changed at next login)")
    except ValidationFailedError as e:
        print(f"Error: {e.detail}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
