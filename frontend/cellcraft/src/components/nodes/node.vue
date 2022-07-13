<template>
  <div>
      <select df-select_type @change="onMyChange">
            <option value="">Select</option>
            <option value="1">Select type 1</option>
            <option value="2">Select type 2</option>
            <option value="3">Select type 3</option>
          </select>
          <div v-if="select_type">
            <p>Select type: {{ select_type }}</p>
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
      title: 'consoleLog'
    }
  },
  props: {
    select_type: String
  },
  methods: {
    onMyChange () {
      console.log(`New select type: ${this.select_type}`)
    }
  }
}
</script>

<style>

</style>
