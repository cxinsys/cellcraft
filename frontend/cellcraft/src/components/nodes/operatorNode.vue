<template>
  <div>
    <nodeHeader v-bind:title="title" />
    <div class="inputBox">
        <input df-operator class="inputBox__choice" type="radio" id="inputBox__choice" name="operator" value="plus">
        <label class="inputBox__text" for="inputBox__choice">plus</label>

        <input df-operator class="inputBox__choice" type="radio" id="inputBox__choice" name="operator" value="minus">
        <label class="inputBox__text" for="inputBox__choice">minus</label>

        <input df-operator class="inputBox__choice" type="radio" id="inputBox__choice" name="operator" value="multiplied">
        <label class="inputBox__text" for="inputBox__choice">multiplied</label>

        <input df-operator class="inputBox__choice" type="radio" id="inputBox__choice" name="operator" value="divided">
        <label class="inputBox__text" for="inputBox__choice">divided</label>
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
      title: 'Operation'
    }
  },
  methods: {

  }
}
</script>

<style>
.inputBox{
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}
.inputBox__choice{
  width:15px;
  height:15px;
  margin: 5px 0;
}
.inputBox__text{
  font-size: 15px;
}
</style>
