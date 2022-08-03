<template>
  <div>
    <!-- <nodeHeader v-bind:title="title" /> -->
    <div class="buttonBox">
        Data Table
    </div>
  </div>
</template>

<script>
// import nodeHeader from '@/components/nodes/nodeHeader.vue'

export default {
  // components: {
  //   nodeHeader
  // },
  data () {
    return {
      title: 'Data Table'
    }
  },
  methods: {

  }
}
</script>

<style>
.buttonBox{
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
  font-size: 1rem;
  font-weight: bold;
}
</style>
