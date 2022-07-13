<template>
  <div>
    <nodeHeader v-bind:title="title" />
    <div>
        <input type="number" />
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
      title: 'Input Num'
    }
  },
  methods: {

  }
}
</script>

<style>

</style>
