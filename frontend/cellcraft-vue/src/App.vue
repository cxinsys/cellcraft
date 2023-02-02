<template>
  <div
    class="app"
    v-bind:class="{
      loginPage__layout: loginPage,
      workflowPage__layout: workflowPage,
    }"
  >
    <header
      class="header-component"
      v-bind:class="{
        loginPage__header: loginPage,
        workflowPage__header: workflowPage,
      }"
    >
      <CellcraftHeader></CellcraftHeader>
    </header>
    <main
      class="main-component"
      v-bind:class="{
        loginPage__main: loginPage,
        workflowPage__main: workflowPage,
      }"
    >
      <router-view></router-view>
    </main>
    <footer
      class="footer-component"
      v-bind:class="{
        loginPage__footer: loginPage,
        workflowPage__footer: workflowPage,
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
    };
  },
  watch: {
    $route(to) {
      console.log(to.path);
      if (to.path === "/main") {
        this.mainPage = true;
        this.loginPage = false;
        this.workflowPage = false;
      } else if (to.path === "/login" || to.path === "/signup") {
        this.mainPage = false;
        this.loginPage = true;
        this.workflowPage = false;
      } else if (to.path === "/workflow") {
        this.mainPage = false;
        this.loginPage = false;
        this.workflowPage = true;
      }
    },
  },
};
</script>

<style>
.app {
  width: 100vw;
  height: 191vh;
}
.loginPage__layout {
  height: 100vh;
}
.workflowPage__layout {
  height: 100vh;
}
.header-component {
  width: 100%;
  /* height: 4%; */
  height: 3%;
}
.loginPage__header {
  display: none;
}
.workflowPage__header {
  /* height: 8.525%; */
  height: 4.6%;
}
.main-component {
  width: 100%;
  /* height: 89%; */
  height: 91.5%;
}
.loginPage__main {
  width: 100%;
  height: 100%;
}
.workflowPage__main {
  /* height: 91.475%; */
  height: 95.4%;
}
.footer-component {
  width: 100%;
  height: 5.5%;
  background: #fff;
  border-top: 1px solid #e1e1e1;
}
.loginPage__footer {
  display: none;
}
.workflowPage__footer {
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
</style>
