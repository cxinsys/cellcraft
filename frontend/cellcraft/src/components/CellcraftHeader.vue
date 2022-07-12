<template>
  <div>
    <section class = "sidebar" id="sidebar" v-bind:class="{open: S_isActive}">

      <div class="logo-details" style="display: flex;align-items: center; justify-content: space-between;">
        <div class="logo_name" style="margin-left: 1.5vw;">
          <a href="" @click="redirect">LOGO</a>
        </div>
        <i class="fa-solid fa-bars" style="margin-right: .5vw;" id="btn" @click="openSidebar"></i>
      </div>

      <ul class="nav-list">
        <li>
          <div class="mainMenu" @click="openMenu1">
            <i class="fa fa-home"></i>
            <span class="links_name">Workflows</span>
          </div>
          <ul class="sub-menu" v-bind:class="{open: M1_isActive}">
            <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
          </ul>
        </li>

        <li>
          <div class="mainMenu" @click="openMenu2">
            <i class="fa fa-user"></i>
            <span class="links_name">Files</span>
          </div>
          <ul class="sub-menu" v-bind:class="{open: M2_isActive}">
            <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
          </ul>
        </li>

        <li>
          <div class="mainMenu" @click="openMenu3">
            <i class="fa fa-envelope-open icon"></i>
            <span class="links_name">Data Sets</span>
          </div>
          <ul class="sub-menu" v-bind:class="{open: M3_isActive}">
            <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
          </ul>
        </li>

        <li>
          <div class="mainMenu" @click="openMenu3">
            <i class="fa fa-envelope-open icon"></i>
            <span class="links_name">Plots</span>
          </div>
          <ul class="sub-menu" v-bind:class="{open: M3_isActive}">
            <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
          </ul>
        </li>

      </ul>
    </section>

    <section class="home-section" id="title">
      <section class="section1" style="margin-left: 1vw;">
        <button @click="openSidebar">
          <i class="fa fa-solid fa-bars" id="btn"></i>
        </button>

        <button>
          <div id = "logo">
            <a href="" @click="redirect">LOGO</a>
          </div>
        </button>
      </section>

      <section class="section2" style="display: flex;">

        <!-- <li class="menu">
          <a href="">
            <i class="fa fa-home"></i>
            <span class="links_name">Home</span>
          </a>
          <ul class="sub-menu" >
            <li><i class="fa-solid fa-heart"></i><a href="">Content1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Content2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Content3</a></li>
          </ul>
        </li>

        <li class="menu">
          <a href="">
            <i class="fa fa-user"></i>
            <span class="links_name">Content</span>
          </a>
          <ul class="sub-menu">
            <li><i class="fa-solid fa-heart"></i><a href="">Content1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Content2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Content3</a></li>
          </ul>
        </li>

        <li class="menu">
          <a href="">
            <i class="fa fa-envelope-open"></i>
            <span class="links_name">Contact</span>
          </a>
          <ul class="sub-menu">
            <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
          </ul>
        </li> -->

        <template v-if="isUserLogin">
          <li class="menu">
            <a href="#" @click="logoutUser" class="links_name">logout</a>
          </li>
          <li class="menu">
            <a href="#" @click="getCurrentUser" class="links_name">UserInfo</a>
          </li>
        </template>

        <template v-else>
          <li class="menu">
            <router-link to="/login" class="links_name">Login</router-link>
          </li>
        </template>

      </section>
    </section>
  </div>
</template>

<script>
import { getUser } from '@/api/index'
import { deleteCookie } from '@/utils/cookies'

export default {
  data () {
    return {
      S_isActive: false, // sidebar
      M1_isActive: false, // subMenu
      M2_isActive: false,
      M3_isActive: false
    }
  },
  computed: {
    isUserLogin () {
      return this.$store.getters.isLogin
    }
  },
  methods: {
    redirect () {
      this.$router.push('/')
    },
    logoutUser () {
      this.$store.commit('clearToken')
      deleteCookie()
      this.$router.push('/')
    },
    async getCurrentUser () {
      const userInfo = await getUser()
      this.$store.commit('setUserInfo', userInfo.data)
      console.log(userInfo.data)
      return userInfo.data
    },
    openSidebar () {
      this.S_isActive = !this.S_isActive
      this.M1_isActive = false
      this.M2_isActive = false
      this.M3_isActive = false
    },
    openMenu1 () {
      this.M1_isActive = !this.M1_isActive
    },
    openMenu2 () {
      this.M2_isActive = !this.M2_isActive
    },
    openMenu3 () {
      this.M3_isActive = !this.M3_isActive
    }
  }

}
</script>

<style>
  @import '../../css/style.css';
  #title {
    background-color: skyblue;
    display: flex;
    justify-content: space-between;
    align-items: center;
    height: 6vh;
  }
  #logo {
    width: 2vw;
    height: auto;
  }
  button {
    border: 0;
    background-color: transparent;
  }

  .sidebar {
    position: fixed;
    height: 100%;
    background: white;
    z-index: 99;
    transition: all 0.5s ease;
    /* display: none; */
    width: 0;
    box-sizing: border-box;
  }

  .sidebar .nav-list {
    display: none;
    height: 6vh;
    list-style: none;
  }
  .sidebar.open .nav-list {
    display: block;
  }

  .sidebar.open{
    transition: all 0.5s ease;
    width: 15vw;
    box-sizing: border-box;
  }

  .sidebar .logo-details{
    height: 6vh;
    position: relative;
    background-color: rgb(123, 123, 123);
    color: white;
  }
  .sidebar .logo-details .logo_name{
    color: white;
    font-size: 1.5vw;
    font-weight: 600;
    opacity: 0;
    transition: all 0.5s ease;
    display: none;
  }
  .sidebar.open .logo-details .logo_name{
    display: block;
    opacity: 1;
  }
  .sidebar .logo-details #btn{
    position: absolute;
    top: 50%;
    right: 0;
    transform: translateY(-50%);
    transition: all 0.4s ease;
    font-size: 1.3vw;
    text-align: center;
    cursor: pointer;
    transition: all 0.5s ease;
  }
  .sidebar.open .logo-details #btn{
    text-align: right;
  }

  .sidebar li{
    position: relative;
    display: list-item;
  }

  .sidebar li div{
    height: 100%;
    width: 90%;
    align-items: center;
    text-decoration: none;
    transition: all 0.4s ease;
    padding: 1.5vh 0;
    padding-left: 1vw;
    padding-right: .5vw;
    /* border-bottom: rgb(206, 206, 206) .1vw solid; */
    border-top: rgb(206, 206, 206) .1vw solid;
  }
  .sidebar li div:hover{
    background-color: rgb(206, 206, 206);
    color: black;
  }

  .sidebar li ul {
    position: relative;
    display: none;
    /* margin-top: 2vh; */
  }

  .sidebar.open .sub-menu.open{
    display: block;
  }

  .sidebar.open .sub-menu.open li{
    padding: 1vw;
    /* padding-left: 2vw; */
    margin-left: 1vw;
    width: 70%;
  }

  .sidebar.open .sub-menu.open li a{
    text-decoration: none;
    color: black;
    padding-left: 1vw;
  }

  .sidebar.open .sub-menu.open li:hover{
    background-color: rgb(206, 206, 206);
    color: black;
  }

  .sidebar li a .links_name, .sidebar li a .fa{
    color: black;
    font-size: 1.2vw;
    font-weight: 400;
    white-space: nowrap;
    opacity: 0;
    pointer-events: none;
    transition: 0.4s;
    margin-bottom: 1vh;
  }

  .sidebar.open li a .links_name, .sidebar li a .fa{
    opacity: 1;
    pointer-events: auto;
  }

  .sidebar li a:hover .links_name,
  .sidebar li a:hover i{
    transition: all 0.5s ease;
    color: #473cbf;
  }

  .home-section{
    position: relative;
    background: #E4E9F7;
    transition: all 0.5s ease;
    z-index: 2;
  }

  .sidebar.open ~ .home-section .section1{
    opacity: 0;
  }

  .section2 .menu > a{
    text-decoration-line: none;
    color: white;
  }
  .section2 a{
    text-decoration-line: none;
    color: white;
  }

  .section2 .menu{
    list-style: none;
    position: relative;
    margin-bottom: -18vh;
  }

  .section2 .links_name {
    margin-right: 1vw;
  }

  .section2 .menu>a:hover {
    text-decoration: underline;
    text-underline-offset: .3vw;
  }

  .section2 .menu:hover ul{
    opacity: 1;
  }

  .section2 li ul{
    position: relative;
    opacity: 0;
    margin-top: 2vh;
    border: white solid;
    background-color: white;
    border-radius: .6vw;
    padding: 1vw;
    list-style: none;
  }

  .section2 li ul li {
    margin-bottom: 2vh;

  }

  .section2 li ul li a{
    color: black;

  }

  .section2 li ul li i {
    margin-right: .5vw;
  }

</style>
