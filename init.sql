-- 데이터베이스 사용자 생성 (idempotent - 이미 존재하면 무시)
DO $$
BEGIN
   IF NOT EXISTS (SELECT FROM pg_catalog.pg_roles WHERE rolname = 'cellcraft_admin') THEN
      CREATE USER cellcraft_admin WITH PASSWORD 'cellcraft_admin';
   END IF;
END
$$;

-- 데이터베이스 생성 (이미 존재하면 오류 무시)
-- Note: CREATE DATABASE cannot be executed from a function or multi-command script
SELECT 'CREATE DATABASE cellcraft OWNER cellcraft_admin'
WHERE NOT EXISTS (SELECT FROM pg_database WHERE datname = 'cellcraft')\gexec

-- public 스키마 권한 부여
GRANT USAGE ON SCHEMA public TO cellcraft_admin;
GRANT CREATE ON SCHEMA public TO cellcraft_admin;

-- 데이터베이스 권한 부여
GRANT ALL PRIVILEGES ON DATABASE cellcraft TO cellcraft_admin;
GRANT SELECT, INSERT, UPDATE, DELETE ON ALL TABLES IN SCHEMA public TO cellcraft_admin;