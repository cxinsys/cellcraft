<style scoped>
  .files__header{
    width: 100vw;
    height: 5vh;
    display: flex;
    align-items: center;
    margin-top: 50px;
    border-bottom: 1px #a1a1a1 solid;
  }
  .l-col{
    width: 25%;
    height: 100%;

    display: flex;
    align-items: center;
    justify-content: center;
  }
  .l-col:first-child{
    width: 50%;
    height: 100%;

    display: flex;
    align-items: center;
    justify-content: center;
  }
  .files__header__img{
    margin-left: 5px;
    width: 15px;
  }
  .font-setting{
    font-size: 1vw;
    font-weight: bold;
  }
  .files__contents{
    width: 100vw;
    height: 75vh;

    margin-top: 20px;
  }
  .files__contents__item{
    display: flex;
    align-items: center;

    width: 100%;
    height: 40px;
  }
</style>

<template>
  <section>
    <div class="files__header">
      <div class="l-col">
        <div class="font-setting"> Name </div>
        <img class="files__header__img" src="@/assets/descending-filter-desc.png" alt="">
      </div>
      <div class="l-col">
        <div class="font-setting"> Size </div>
      </div>
      <div class="l-col">
        <div class="font-setting"> Date </div>
      </div>
    </div>
    <ul class="files__contents">
      <li class="files__contents__item" v-for="(file, idx) in files_list" :key="idx">
        <div class="l-col">{{file.file_name}}</div>
        <div class="l-col">{{file.file_size}}</div>
        <div class="l-col">{{file.created_at}}</div>
      </li>
    </ul>
  </section>
</template>

<script>
import { getFiles } from '@/api/index'

export default {
  data () {
    return {
      files_list: []
    }
  },
  async mounted () {
    const fileList = await getFiles()
    console.log(fileList.data)
    this.files_list = fileList.data
  }
}
</script>
