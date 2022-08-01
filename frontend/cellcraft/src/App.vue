<template>
  <div id="app">
    <CellcraftHeader></CellcraftHeader>
    <div class="app-contents">
      <router-view></router-view>
    </div>
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
