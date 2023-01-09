<template>
  <div class="app" v-bind:class="{loginPage__layout: loginPage}">
    <header class="header-component" v-bind:class="{loginPage__header: loginPage}">
      <CellcraftHeader></CellcraftHeader>
    </header>
    <main class="main-component" v-bind:class="{loginPage__main: loginPage}">
      <router-view></router-view>
    </main>
    <footer class="footer-component" v-bind:class="{loginPage__footer: loginPage}">
      <div class="copyright__txt">
        CELLCRAFT © 2023. All rights reserved
      </div>
    </footer>
  </div>

</template>

<script>
import CellcraftHeader from '@/components/CellcraftHeader.vue'

export default {
  components: {
    'CellcraftHeader': CellcraftHeader
  },
  data () {
    return {
      mainPage: true,
      loginPage: false,
      workflowPage: false
    }
  },
  watch: {
    '$route' (to, from) {
      console.log(to.path)
      if (to.path === '/main') {
        this.mainPage = true
        this.loginPage = false
        this.workflowPage = false
      } else if (to.path === '/login') {
        this.mainPage = false
        this.loginPage = true
        this.workflowPage = false
      } else if (to.path === '/workflow') {
        this.mainPage = false
        this.loginPage = false
        this.workflowPage = true
      }
    }
  }
}
</script>

<style>
.app{
  width: 100vw;
  height: 155vh;
}
.loginPage__layout{
  height: 100vh;
}
.workflowPage__layout{
  height: 100vh;
}
.header-component{
  width: 100%;
  height: 5.5%;
}
.loginPage__header{
  display: none;
}
.workflowPage__header{
}
.main-component{
  width: 100%;
  height: 89%;
}
.loginPage__main{
  width: 100%;
  height: 100%;
}
.workflowPage__main{
}
.footer-component{
  width: 100%;
  height: 5.5%;
  background: #fff;
  border-top: 1px solid #e1e1e1;
}
.loginPage__footer{
  display: none;
}
.workflowPage__footer{
}
.copyright__txt{
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;

  font-family: 'NanumGothic';
  font-style: normal;
  font-weight: 400;
  font-size: 16px;
  line-height: 16px;
  /* identical to box height */

  color: #5E5E5E;
}
</style>
