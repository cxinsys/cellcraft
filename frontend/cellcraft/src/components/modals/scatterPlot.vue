<template>
<div>
    <h1>Plot</h1>
    <img id="plot" src=''>
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
      node_name: 'ScatterPlot',
      img: null
    }
  },
  async mounted () {
    const filename = {'filename': `${this.node_name}_${this.file_name}`}
    console.log(filename)
    const PlotResult = await getResult(filename)
    console.log(PlotResult)
    this.img = PlotResult.data
    document.getElementById('plot').src = 'data:image/png;base64,' + PlotResult.data
  }
}
</script>

<style>

</style>
