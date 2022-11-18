<template>
  <div class="layout">
    <div>{{ data }}</div>
  </div>
</template>

<script>
import { getResult } from '@/api/index'
export default {
  props: {
    file_name: null
  },
  data () {
    return {
      node_name: 'Plot',
      data: []
    }
  },
  async mounted () {
    const filename = { filename: `${this.node_name}_${this.file_name}` }
    console.log(filename)
    const PlotResult = await getResult(filename)
    console.log(PlotResult)
    this.data = PlotResult.data
    // document.getElementById('plot').src =
    //   'data:image/png;base64,' + PlotResult.data
  }
}
</script>

<style>
#plot {
  width: 100%;
  height: 100%;
  object-fit: contain;
}
.layout {
  width: 100%;
  height: 100%;
}
</style>
