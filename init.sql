-- 데이터베이스 사용자 생성
-- CREATE USER cellcraft WITH PASSWORD 'cellcraft';

-- 데이터베이스 생성
CREATE DATABASE cellcraft OWNER cellcraft;

-- public 스키마 권한 부여
GRANT USAGE ON SCHEMA public TO cellcraft;
GRANT CREATE ON SCHEMA public TO cellcraft;

-- 데이터베이스 권한 부여
GRANT ALL PRIVILEGES ON DATABASE cellcraft TO cellcraft;
GRANT SELECT, INSERT, UPDATE, DELETE ON ALL TABLES IN SCHEMA public TO cellcraft;