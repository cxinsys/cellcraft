<template>
  <div class="layout">
    <nodeHeader v-bind:title="title" />
    <div class="inputBox">
        <input df-input_url class="inputBox__input" type="url" placeholder="Request URL" />
    </div>
    <div class="selectBox">
      <select class="selectBox__select" df-select_type v-model="select_type">
      <option value="">Select</option>
      <option value="GET">GET</option>
      <option value="POST">POST</option>
    </select>
    </div>
    <div class="infoReq" v-if="select_type">
      <p>HTTP Req: {{ select_type }}</p>
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
      title: 'Select Type',
      select_type: ''
    }
  },
  methods: {

  }
}
</script>

<style>
.layout{
  width:100%;
  height:100%;
}
.inputBox ,.selectBox{
  display: flex;
  justify-content: center;
  align-items: center;
  margin: 10px;
}
.selectBox__select{
  width: 100%;
  height: auto;
}
.infoReq{
  display: flex;
  justify-content: start;
  align-items: center;
  margin: 20px;
  width: 150px;
  height: auto;
}
</style>
