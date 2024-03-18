<template>
  <div
    class="app"
    v-bind:class="{
      main__app: mainPage,
    }"
  >
    <header
      class="header-component"
      v-bind:class="{
        loginPage__header: loginPage,
        workflowPage__header: workflowPage,
        filesPage__header: filesPage,
        datasetsPage__header: datasetsPage,
        adminPage__header: adminPage,
        tutorialPage__header: tutorialPage,
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
        datasetsPage__main: datasetsPage,
        adminPage__main: adminPage,
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
        datasetsPage__footer: datasetsPage,
        adminPage__footer: adminPage,
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
      datasetsPage: false,
      adminPage: false,
      tutorialPage: false,
      pluginsPage: false,
    };
  },
  watch: {
    $route(to, from) {
      if (from.path.includes("/workflow") && !to.path.includes("/workflow")) {
        this.$store.commit("clearWorkflow");
        this.$store.commit("clearNodes");
        this.$store.commit("clearLinkedNodes");
        this.$store.commit("clearTitle");
        this.$store.commit("clearThumbnail");
        this.$store.commit("clearCurrentNode");
      }

      if (from.path.includes("/login")) {
        this.$store.commit("clearWorkflow");
        this.$store.commit("clearNodes");
        this.$store.commit("clearLinkedNodes");
        this.$store.commit("clearTitle");
        this.$store.commit("clearThumbnail");
        this.$store.commit("clearCurrentNode");
      }

      if (to.path.includes("/main")) {
        this.mainPage = true;
        this.loginPage = false;
        this.workflowPage = false;
        this.filesPage = false;
        this.datasetsPage = false;
        this.adminPage = false;
        this.tutorialPage = false;
      } else if (to.path.includes("/login") || to.path.includes("/signup")) {
        this.mainPage = false;
        this.loginPage = true;
        this.workflowPage = false;
        this.filesPage = false;
        this.datasetsPage = false;
        this.adminPage = false;
        this.tutorialPage = false;
      } else if (to.path.includes("/workflow")) {
        this.mainPage = false;
        this.loginPage = false;
        this.workflowPage = true;
        this.filesPage = false;
        this.datasetsPage = false;
        this.pluginsPage = false;
        this.adminPage = false;
        this.tutorialPage = false;
      } else if (to.path.includes("/files") || to.path.includes("/projects")) {
        this.mainPage = false;
        this.loginPage = false;
        this.workflowPage = false;
        this.filesPage = true;
        this.datasetsPage = false;
        this.pluginsPage = false;
        this.adminPage = false;
        this.tutorialPage = false;
      } else if (
        to.path.includes("/datasets") ||
        to.path.includes("/plugins")
      ) {
        this.mainPage = false;
        this.loginPage = false;
        this.workflowPage = false;
        this.filesPage = false;
        this.datasetsPage = true;
        this.pluginsPage = true;
        this.adminPage = false;
        this.tutorialPage = false;
      } else if (to.path.includes("/admin")) {
        this.mainPage = false;
        this.loginPage = false;
        this.workflowPage = false;
        this.filesPage = false;
        this.datasetsPage = false;
        this.pluginsPage = false;
        this.adminPage = true;
        this.tutorialPage = false;
      } else if (to.path.includes("/tutorial")) {
        this.mainPage = false;
        this.loginPage = false;
        this.workflowPage = false;
        this.filesPage = false;
        this.datasetsPage = false;
        this.pluginsPage = false;
        this.adminPage = false;
        this.tutorialPage = true;
      }
    },
  },
};
</script>

<style>
.app {
  width: 100vw;
  /* Changed from 100vw to 80vw */
  overflow-x: hidden;
}

.main__app {
  width: calc(100vw - 15px);
  /* Changed from 100vw to 80vw */
}

.header-component {
  width: 100%;
  height: 48px;
  /* Changed from 60px to 48px */
  position: fixed;
  z-index: 999;
  background-color: rgba(0, 0, 0, 0.9);
}

.workflowPage__header,
.filesPage__header,
.datasetsPage__header,
.adminPage__header {
  width: 100%;
  height: 48px;
  /* Changed from 60px to 48px */
  position: relative;
  background-color: rgba(0, 0, 0, 0.9);
}

.tutorialPage__header {
  display: none;
}

.main-component {
  width: 100%;
}

.loginPage__main {
  width: 100%;
  height: 80%;
  /* Assuming you want to change the height proportionally */
  background: #f0f0f0;
}

.workflowPage__main,
.filesPage__main,
.datasetsPage__main,
.adminPage__main {
  height: calc(100vh - 48px);
  /* Changed from 100vh to 80vh and adjusted the subtraction */
}

.footer-component {
  width: 100%;
  height: 80px;
  /* Changed from 80px to 64px */
  background: #fff;
  border-top: 1px solid #e1e1e1;
}

.loginPage__footer,
.workflowPage__footer,
.filesPage__footer,
.datasetsPage__footer {
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

/* You can apply similar calculations for any other properties that need to be scaled down to 80% */
</style>
