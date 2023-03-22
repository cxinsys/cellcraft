<template>
  <div class="app">
    <header
      class="header-component"
      v-bind:class="{
        loginPage__header: loginPage,
        workflowPage__header: workflowPage,
        filesPage__header: filesPage,
      }"
    >
      <CellcraftHeader></CellcraftHeader>
    </header>
    <main
      class="main-component"
      v-bind:class="{
        loginPage__main: loginPage,
        workflowPage__main: workflowPage,
        filesPage__main: filesPage,
      }"
    >
      <router-view></router-view>
    </main>
    <footer
      class="footer-component"
      v-bind:class="{
        loginPage__footer: loginPage,
        workflowPage__footer: workflowPage,
        filesPage__footer: filesPage,
      }"
    >
      <div class="copyright__txt">CELLCRAFT © 2023. All rights reserved</div>
    </footer>
  </div>
</template>

<script>
import CellcraftHeader from "@/components/CellcraftHeader.vue";

export default {
  components: {
    CellcraftHeader: CellcraftHeader,
  },
  data() {
    return {
      mainPage: true,
      loginPage: false,
      workflowPage: false,
      filesPage: false,
    };
  },
  watch: {
    $route(to, from) {
      console.log(to.path, from.path);
      if (from.path === "/workflow") {
        this.$store.commit("clearNodes");
        this.$store.commit("clearLinkedNodes");
      }

      if (to.path === "/main") {
        this.mainPage = true;
        this.loginPage = false;
        this.workflowPage = false;
        this.filesPage = false;
      } else if (to.path === "/login" || to.path === "/signup") {
        this.mainPage = false;
        this.loginPage = true;
        this.workflowPage = false;
        this.filesPage = false;
      } else if (to.path === "/workflow") {
        this.mainPage = false;
        this.loginPage = false;
        this.workflowPage = true;
        this.filesPage = false;
      } else if (to.path === "/files" || to.path === "/projects") {
        this.mainPage = false;
        this.loginPage = false;
        this.workflowPage = false;
        this.filesPage = true;
      }
    },
  },
};
</script>

<style>
.app {
  width: 100vw;
  height: 100vh;
  overflow-x: hidden;
  /* height: 191vh; */
}

.header-component {
  width: 100%;
  /* height: 3.5rem; */
  height: 60px;
  position: fixed;
  z-index: 999;
  background-color: rgba(0, 0, 0, 0.7);
}

/* .header-component, */
.workflowPage__header,
.filesPage__header {
  width: 100%;
  height: 60px;
  position: relative;
  background-color: rgba(0, 0, 0, 0.9);
}

.loginPage__header {
  display: none;
}
.filesPage__header {
  /* border-bottom: 1px solid #e1e1e1; */
}
.main-component {
  width: 100%;
}
.loginPage__main {
  width: 100%;
  height: 100%;
}
.workflowPage__main {
  /* border-top: 1px solid #e1e1e1; */
}
.workflowPage__main,
.filesPage__main {
  height: calc(100vh - 56px);
}
.footer-component {
  width: 100%;
  height: 80px;
  background: #fff;
  border-top: 1px solid #e1e1e1;
}
.loginPage__footer,
.workflowPage__footer,
.filesPage__footer {
  display: none;
}
.copyright__txt {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;

  font-family: "NanumGothic";
  font-style: normal;
  font-weight: 400;
  font-size: 16px;
  line-height: 16px;
  /* identical to box height */

  color: #5e5e5e;
}

@media (prefers-color-scheme: dark) {
  /* .header-component {
    background-color: rgba(0, 0, 0, 0.7);
  } */
  /* .workflowPage__main {
    border-top: 1px solid #404040;
  } */
}
</style>
