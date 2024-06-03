const { defineConfig } = require("@vue/cli-service");

module.exports = defineConfig({
  transpileDependencies: true,
  runtimeCompiler: true,
  css: {
    extract: true, // CSS를 별도의 파일로 추출
    sourceMap: false, // 소스 맵 비활성화
    loaderOptions: {
      css: {},
      postcss: {},
    },
  },
  configureWebpack: {
    optimization: {
      splitChunks: {
        chunks: 'all',
      },
    },
  },
});
