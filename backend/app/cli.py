"""Operator CLI for on-premise installs (#138).

Run inside the backend container, from the backend root:

    docker compose -f <compose-file> exec backend \
        /opt/micromamba/envs/snakemake/bin/python -m app.cli reset-password --email <email>

Commands:
    create-admin    create an admin account (temporary password unless --prompt-password is given)
    reset-password  issue a temporary password; the user must change it at next login
                    (inactive accounts stay inactive unless --activate is given)
"""
import argparse
import getpass
import logging
import sys

from app.core.bootstrap import (
    DEFAULT_ADMIN_USERNAME, _validate_email, generate_password, remove_initial_password_file,
    validate_new_password,
)
from app.core.exceptions import ValidationFailedError
from app.core.security import get_password_hash
from app.db.base_registry import Base  # noqa: F401  (registers every model for relationships)
from app.db.session import get_db_session
from app.user import models


def create_admin(db, *, email: str, username: str, password: str = None) -> str:
    email = _validate_email(email)
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


def reset_password(db, *, email: str, activate: bool = False) -> str:
    user = db.query(models.User).filter_by(email=email).first()
    if not user:
        raise ValidationFailedError(f"User not found: {email}")
    password = generate_password()
    user.hashed_password = get_password_hash(password)
    user.must_change_password = True
    if activate:
        user.is_active = True
    db.commit()
    # 이전 초기 비밀번호 안내 파일은 더 이상 유효하지 않다
    remove_initial_password_file(email)
    return password


def _quiet_sql_logging() -> None:
    """운영자 출력이 SQL echo 로그에 묻히지 않도록 한다."""
    from app.db.session import engine
    engine.echo = False
    logging.getLogger("sqlalchemy.engine").setLevel(logging.WARNING)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m app.cli", description="CellCraft operator commands")
    commands = parser.add_subparsers(dest="command", required=True)

    create = commands.add_parser("create-admin", help="create an admin account")
    create.add_argument("--email", required=True)
    create.add_argument("--username", default=DEFAULT_ADMIN_USERNAME)
    create.add_argument("--prompt-password", action="store_true",
                        help="enter the password interactively instead of generating a temporary one")

    reset = commands.add_parser("reset-password", help="issue a temporary password")
    reset.add_argument("--email", required=True)
    reset.add_argument("--activate", action="store_true", help="also re-activate a deactivated account")

    args = parser.parse_args(argv)
    _quiet_sql_logging()
    try:
        password = None
        if args.command == "create-admin" and args.prompt_password:
            password = getpass.getpass("Password: ")
            if password != getpass.getpass("Confirm password: "):
                raise ValidationFailedError("Passwords do not match")
        with get_db_session() as db:
            if args.command == "create-admin":
                temporary = create_admin(db, email=args.email, username=args.username, password=password)
                print(f"Admin account created: {args.email}")
                if password is None:
                    print(f"Temporary password: {temporary}  (must be changed at first login)")
            else:
                temporary = reset_password(db, email=args.email, activate=args.activate)
                print(f"Temporary password for {args.email}: {temporary}  (must be changed at next login)")
    except ValidationFailedError as e:
        print(f"Error: {e.detail}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
