<template>
  <div>
    <div class="inputBox">
        file
    </div>
  </div>
</template>

<script>
export default {
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
.inputBox{
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
  font-size: 1rem;
  font-weight: bold;

}
.inputBox__input{
  width: 100px;
}
</style>
