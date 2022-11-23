<template>
  <div class="layout">
    <div class="options">
      <div class="options__column">
        <div class="options__name"> Values: </div>
        <!-- <div class="options__button" @click="optionClick"> {{ this.value }}</div> -->
        <input type="text" class="options__button" name="value" v-model="value" @click="optionClick">
        <ul class="options__menu">
          <li class="options__content" v-for="(header, idx) in headers" :key="idx" @click="menuClick">
            {{header}}
          </li>
        </ul>
      </div>
      <div class="options__column">
        <div class="options__name"> Annotations: </div>
        <!-- <div class="options__button" @click="optionClick"> {{ this.annotation }}</div> -->
        <input type="text" class="options__button" name="annotation" v-model="annotation" @click="optionClick">
        <ul class="options__menu">
          <li class="options__content" v-for="(header, idx) in headers" :key="idx" @click="menuClick">
            {{header}}
          </li>
        </ul>
      </div>
    </div>
  </div>
</template>

<script>
import { getResult } from '@/api/index'
export default {
  props: {
    file_name: null
  },
  data () {
    return {
      node_name: 'Plot',
      data: [],
      headers: [],
      value: '',
      annotation: ''
    }
  },
  async mounted () {
    const filename = { filename: `${this.node_name}_${this.file_name}` }
    // console.log(filename)
    const PlotResult = await getResult(filename)
    // console.log(PlotResult.data)
    this.data = PlotResult.data
    this.headers = Object.keys(PlotResult.data).slice(1)
    this.value = this.headers[0]
    this.annotation = this.headers[1]
  },
  methods: {
    optionClick (event) {
      event.target.style.opacity = '0.6'
      event.target.nextElementSibling.style.display = 'block'
    },
    menuClick (event) {
      console.dir(event.target.parentElement.previousElementSibling)
      const menuContent = event.target.innerText
      event.target.parentElement.style.display = 'none'
      event.target.parentElement.previousElementSibling.style.opacity = '1'
      if (event.target.parentElement.previousElementSibling.previousElementSibling.innerText.includes('Annotation')) {
        this.annotation = menuContent
      } else if (event.target.parentElement.previousElementSibling.previousElementSibling.innerText.includes('Value')) {
        this.value = menuContent
      }
    }
  }
}
</script>

<style>
.options, .options__column, .options__name, .options__button, .options__menu, .options__content{
  box-sizing: border-box;
  border-radius: 5px;
}
.layout {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  position: relative;
}
.options{
  width: 40%;
  height: 100%;
  padding: 10px;
  display: flex;
  flex-direction: column;
  align-items: center;
  background: rgb(235, 235, 235);
  position: absolute;
  left: 0;
}
.options__column{
  width: 100%;
  height: 5%;
  margin-bottom: 15px;
  display: flex;
  align-items: center;
  position: relative;
  background: rgb(210, 210, 210);
}
.options__name{
 width: 35%;
 height: 100%;
 padding: 5px;
 display: flex;
 align-items: center;
}
.options__button, .options__menu{
  width: 65%;
  height: 100%;
  padding: 5px;
  display: flex;
  align-items: center;
  justify-content: flex-start;
  background: white;
}
.options__menu{
  display: none;
  height: auto;
  flex-direction: column;
  position: absolute;
  box-sizing: border-box;
  border-radius: 7px;
  top: 100%;
  right: 0;
  z-index: 1;
}
.options__content{
  margin-bottom: 4px;
}
.options__content:hover{
  background: rgb(125, 125, 255);
}
</style>
