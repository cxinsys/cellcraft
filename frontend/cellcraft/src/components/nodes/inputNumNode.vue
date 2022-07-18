<template>
  <div>
    <nodeHeader v-bind:title="title" />
    <div class="inputBox">
        <input df-input_num class="inputBox__input" type="number" value="0" />
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
.inputBox__input{
  width: 100px;
}
</style>
