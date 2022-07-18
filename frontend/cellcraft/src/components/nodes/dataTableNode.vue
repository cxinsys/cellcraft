<template>
  <div>
    <nodeHeader v-bind:title="title" />
    <div class="buttonBox">
        <button class="buttonBox__button">Data Table</button>
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
      title: 'Data Table'
    }
  },
  methods: {

  }
}
</script>

<style>
.buttonBox{
  display: flex;
  align-items: center;
  justify-content: center;
}
.buttonBox__button{
  display: flex;
  align-items: center;
  justify-content: center;
  width: 100px;
  height: 50px;
  border-radius: 8px;
  border: 2px solid #494949;
  line-height:40px;
  padding: 10px;
  margin: 10px 0px;
  cursor: pointer;
}
</style>
