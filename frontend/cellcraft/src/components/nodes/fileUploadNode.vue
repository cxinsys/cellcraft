<template>
  <div>
    <nodeHeader v-bind:title="title" />
    <div class="inputBox">
        <input df-file class="inputBox__input" type="file" />
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
      title: 'File Upload'
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
