<template>
  <div id="app">
    <header class="header-component">
      <CellcraftHeader></CellcraftHeader>
    </header>
    <main class="main-component">
      <router-view></router-view>
    </main>
    <footer class="footer-component">
      <div class="copyright__txt">
        Cell Craft © 2023. All rights reserved
      </div>
    </footer>
  </div>

</template>

<script>
import CellcraftHeader from '@/components/CellcraftHeader.vue'

export default {
  components: {
    'CellcraftHeader': CellcraftHeader
  }
}
</script>

<style>
#app{
  width: 100vw;
  height: 155vh;
}
.header-component{
  width: 100%;
  height: 5.5%;
}
.main-component{
  width: 100%;
  height: 89%;
}
.footer-component{
  width: 100%;
  height: 5.5%;
  background: #fff;
  border-top: 1px solid #e1e1e1;
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
