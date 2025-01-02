-- 데이터베이스 사용자 생성
CREATE USER testdb WITH PASSWORD 'testdb';

-- 데이터베이스 생성
CREATE DATABASE cellcraft OWNER testdb;

-- public 스키마 권한 부여
GRANT USAGE ON SCHEMA public TO testdb;
GRANT CREATE ON SCHEMA public TO testdb;

-- 데이터베이스 권한 부여
GRANT ALL PRIVILEGES ON DATABASE cellcraft TO testdb;
GRANT SELECT, INSERT, UPDATE, DELETE ON ALL TABLES IN SCHEMA public TO testdb;