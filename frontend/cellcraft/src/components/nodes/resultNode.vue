<template>
  <div>
    <nodeHeader v-bind:title="title" />
    <div>
        <button class="resultButton">Result</button>
    </div>
  </div>
</template>

<script>
import nodeHeader from '@/components/nodes/nodeHeader.vue'

export default {
  components: {
    nodeHeader
  },
  data () {
    return {
      title: 'Result'
    }
  },
  methods: {

  }
}
</script>

<style>
.resultButton{
    width: 80px;
    height: 20px;
}
</style>
