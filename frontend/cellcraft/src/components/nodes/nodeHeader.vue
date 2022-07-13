<template>
  <div class="header">{{ title }}</div>
</template>

<script>
export default {
  props: ['title']
}
</script>

<style>
.header{
  background: #95eaff;
  margin-top: -15px;
  margin-left: -15px;
  margin-right: -15px;
  padding: 10px 15px;
  margin-bottom: 15px;
}
</style>
