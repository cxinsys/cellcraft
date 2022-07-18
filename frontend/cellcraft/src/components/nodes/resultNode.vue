<template>
  <div>
    <nodeHeader v-bind:title="title" />
    <div class="resultBox">
        <button class="resultBox__Button">Result</button>
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
.resultBox{
  display: flex;
  align-items: center;
  justify-content: center;
}
.resultBox__Button{
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
