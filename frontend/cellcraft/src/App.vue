<template>
  <div id="app">
    <CellcraftHeader></CellcraftHeader>
    <div class="app-contents">
      <router-view></router-view>
    </div>
    <footer class="copyright">
      <div class="copyright__txt">
        <!-- Copyright 2022. Cellcraft Co. all rights reserved -->
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
  height: 100vh;
}
.app-contents{
  width: 100%;
  height: 87%;
}
.copyright{
  display: flex;
  align-items: center;
  justify-content: center;
}
</style>
