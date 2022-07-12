<template>
  <div>
    <section class = "sidebar" v-bind:class="{open: S_isActive}">

      <div class="sidebar__logo">
        <div class="sidebar__logo__name">
          <a href="" @click="redirect">LOGO</a>
        </div>
        <i class="fa-solid fa-bars" @click="openSidebar"></i>
      </div>

      <ul class="sidebar__navList">
        <li class="i-sidebar__navList">
          <div class="sidebar__mainMenu" @click="openMenu1">
            <i class="fa fa-home"></i>
            <span class="sidebar__linksName">Workflows</span>
          </div>
          <ul class="sidebar__subMenu" v-bind:class="{open: M1_isActive}">
            <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
          </ul>
        </li>

        <li class="i-sidebar__navList">
          <div class="sidebar__mainMenu" @click="openMenu2">
            <i class="fa fa-user"></i>
            <span class="sidebar__linksName">Files</span>
          </div>
          <ul class="sidebar__subMenu" v-bind:class="{open: M2_isActive}">
            <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
          </ul>
        </li>

        <li class="i-sidebar__navList">
          <div class="sidebar__mainMenu" @click="openMenu3">
            <i class="fa fa-envelope-open icon"></i>
            <span class="sidebar__linksName">Data Sets</span>
          </div>
          <ul class="sidebar__subMenu" v-bind:class="{open: M3_isActive}">
            <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
            <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
            <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
          </ul>
        </li>

        <li class="i-sidebar__navList">
          <div class="sidebar__mainMenu">
            <i class="fa fa-envelope-open icon"></i>
            <span class="sidebar__linksName">Plots</span>
          </div>

        </li>
      </ul>
    </section>

    <section class="home homeMain">
      <section class="home__section1">
        <button @click="openSidebar">
          <i class="fa fa-solid fa-bars"></i>
        </button>

        <button>
          <div class = "home__logo">
            <a href="" @click="redirect">LOGO</a>
          </div>
        </button>
      </section>

      <section class="home__section2" style="display: flex;">

        <template v-if="isUserLogin">
          <li class="home__menu">
            <a href="#" @click="logoutUser" class="home__linksName">logout</a>
          </li>
          <li class="home__menu">
            <a href="#" @click="getCurrentUser" class="home__linksName">UserInfo</a>
          </li>
        </template>

        <template v-else>
          <li class="home__menu">
            <router-link to="/login" class="home__linksName">Login</router-link>
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

  .sidebar.open{
    transition: all 0.5s ease;
    width: 15vw;
    box-sizing: border-box;
  }

  .sidebar__navList{
    display: none;
    height: 6vh;
    list-style: none;
  }
  .sidebar.open .sidebar__navList{
    display: block;
  }

  .sidebar__logo{
    height: 6vh;
    position: relative;
    background-color: rgb(123, 123, 123);
    color: white;
    display: flex;
    align-items: center;
    justify-content: space-between;
  }
  .sidebar__logo__name{
    color: white;
    font-size: 1.5vw;
    font-weight: 600;
    opacity: 0;
    transition: all 0.5s ease;
    display: none;
    margin-left: 1.5vw;
  }
  .sidebar.open .sidebar__logo__name{
    display: block;
    opacity: 1;
  }
  .sidebar .fa-bars{
    position: absolute;
    top: 50%;
    right: 0;
    transform: translateY(-50%);
    transition: all 0.4s ease;
    font-size: 1.3vw;
    text-align: center;
    cursor: pointer;
    transition: all 0.5s ease;
    margin-right: .5vw;
  }
  .sidebar.open .fa-bars{
    text-align: right;
  }

  .i-sidebar__navList{
    position: relative;
    display: list-item;
  }

  .sidebar__mainMenu{
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
  .sidebar__mainMenu:hover{
    background-color: rgb(206, 206, 206);
    color: black;
  }

  .sidebar__subMenu {
    position: relative;
    display: none;
    /* margin-top: 2vh; */
  }

  .sidebar.open .sidebar__subMenu.open{
    display: block;
  }

  .sidebar.open .sidebar__subMenu.open li{
    padding: 1vw;
    /* padding-left: 2vw; */
    margin-left: 1vw;
    width: 70%;
  }

  .sidebar.open .sidebar__subMenu.open li a{
    text-decoration: none;
    color: black;
    padding-left: 1vw;
  }

  .sidebar.open .sidebar__subMenu.open li:hover{
    background-color: rgb(206, 206, 206);
    color: black;
  }

  /* .sidebar li a .links_name, .sidebar li a .fa{
    color: black;
    font-size: 1.2vw;
    font-weight: 400;
    white-space: nowrap;
    opacity: 0;
    pointer-events: none;
    transition: 0.4s;
    margin-bottom: 1vh;
  } */

  /* .sidebar.open li a .links_name, .sidebar li a .fa{
    opacity: 1;
    pointer-events: auto;
  } */
/*
  .sidebar li a:hover .links_name,
  .sidebar li a:hover i{
    transition: all 0.5s ease;
    color: #473cbf;
  } */

  .home{
    position: relative;
    background: #E4E9F7;
    transition: all 0.5s ease;
    /* z-index: 2; */
    background-color: white;
    display: flex;
    justify-content: space-between;
    align-items: center;
    height: 6vh;
  }

  /* .homeMain{
    z-index: 2;
    background-color: white;
    display: flex;
    justify-content: space-between;
    align-items: center;
    height: 6vh;
  } */

  .home__section1{
    margin-left: 1vw;
  }

  .sidebar.open ~ .home__section1{
    opacity: 0;
  }

  .home__menu{
    list-style: none;
    position: relative;
    margin-right: 2vw;
  }

  .home__menu > a{
    text-decoration-line: none;
    color: black;
    border: .2vw black solid;
    padding: 0.5vh 1vw;
  }

  .home__menu>a:hover {
    /* text-decoration: underline;
    text-underline-offset: .3vw; */
    font-weight: bold;
  }

  .home__logo > a{
    text-decoration-line: none;
    color: black;

  }

  /* .section2 .menu:hover ul{
    opacity: 1;
  } */
/*
  .section2 li ul{
    position: relative;
    opacity: 0;
    margin-top: 2vh;
    border: white solid;
    background-color: white;
    border-radius: .6vw;
    padding: 1vw;
    list-style: none;
  } */

  /* .section2 li ul li {
    margin-bottom: 2vh;

  }

  .section2 li ul li a{
    color: black;

  }

  .section2 li ul li i {
    margin-right: .5vw;
  } */

</style>
