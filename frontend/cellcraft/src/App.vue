<template>
  <div id="app">
    <CellcraftSidebar></CellcraftSidebar>
    <CellcraftHeader></CellcraftHeader>
  </div>
</template>

<script>
import CellcraftHeader from './components/CellcraftHeader.vue'
import CellcraftSidebar from './components/CellcraftSidebar.vue'

export default {
  components: {
    'CellcraftHeader': CellcraftHeader,
    'CellcraftSidebar': CellcraftSidebar
  }
}
</script>

<style>

</style>
