const { defineConfig } = require('@vue/cli-service')

module.exports = defineConfig({
  transpileDependencies: true,
  runtimeCompiler: true,
  publicPath: "/",
  devServer: {
    client: {
      overlay: false,
    },
  },
  css: {
    extract: true,
    sourceMap: false,
    loaderOptions: {
      css: {},
      postcss: {},
    },
  },
  configureWebpack: {
    optimization: {
      splitChunks: {
        chunks: "all",
      },
    },
  },
});
