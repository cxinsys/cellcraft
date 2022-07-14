<template>
  <div>
    <nodeHeader v-bind:title="title" />
    <div class="inputBox">
        <input type="radio" id="inputBox__choice" name="plus" value="plus">
        <label for="inputBox__choice">plus</label>

        <input type="radio" id="inputBox__choice" name="minus" value="minus">
        <label for="inputBox__choice">minus</label>

        <input type="radio" id="inputBox__choice" name="multiplied" value="multiplied">
        <label for="inputBox__choice">multiplied</label>

        <input type="radio" id="inputBox__choice" name="divided" value="divided">
        <label for="inputBox__choice">divided</label>
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

</style>
