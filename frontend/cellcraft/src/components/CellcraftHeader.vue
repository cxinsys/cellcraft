<style scoped>

 a{
    text-decoration: none;
    color: black;
  }

  button {
    border: 0;
    background-color: transparent;
  }

  .sidebar {
    position: fixed;
    height: 100%;
    background: white;
    z-index: 5;
    transition: all 0.5s ease;
    /* display: none; */
    width: 0;
    box-sizing: border-box;
    /* position: ; */
  }

  .sidebar.open{
    transition: all 0.5s ease;
    width: 13vw;
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
    background-color: #242F9B;
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
  .sidebar__logoImg{
    height: auto;
    width: 10vw;
    margin-left: -1vw;
    margin-top: .7vh;
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
    height: 2vh;
    width: 100%;
    align-items: center;
    text-decoration: none;
    transition: all 0.4s ease;
    padding: 2.5vh 0vh 2.5vh 1vw;
    border-top: rgb(206, 206, 206) .1vw solid;
    display: flex;
  }

  .mainMenu__icon{
    width: 1.5vw;
    height: auto;
    margin-right: 1vw;
  }

  .sidebar__linksName{
    font-weight: bold;
    text-decoration: none;
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

  .home{
    position: relative;
    background: #E4E9F7;
    transition: all 0.5s ease;
    /* z-index: 2; */
    background-color: white;
    display: flex;
    justify-content: space-between;
    align-items: flex-start;
    height: 55px;
  }

  .home__section1{
    margin-left: 1vw;
    display: flex;
    margin-top: -5px;
  }

  .sidebar.open ~ .home__section1{
    opacity: 0;
  }

  .home__login{
    list-style: none;
    position: relative;
    margin-right: 2vw;
    background-color: black;
    border-radius: .7vw;
    /* border: .2vw black solid; */
    padding: 1vh 0;
    margin-top: 7px;
  }

  /* .home__menu > a{
    text-decoration-line: none;
    color: white;
    background-color: black;
    border-radius: .7vw;
    font-size: 1vw;
    font-weight: bold;
    padding: 0.7vh 1vw;
  } */

  .home__login__text{
    text-decoration-line: none;
    color: white;
    border-radius: .7vw;
    font-size: 1vw;
    font-weight: bold;
    padding: 0.7vh 1vw;
  }

  /* .home__login>a:hover {
    text-decoration: underline;
    text-underline-offset: .3vw;
    font-weight: bold;
  } */

  .sidebar.open ~ .home .home__logo{
    opacity: 0;
  }

  .home__logo > a{
    text-decoration-line: none;
    color: black;

  }
  .home__logoImg {
    height: 6vh;
    margin-top: .7vh;
  }
  /* .home__login{
    color: white;
    text-decoration-line: none;
    font-size: .9vw;
    font-weight: bold;
  } */
  .home__linksName{
    color: white;
    text-decoration-line: none;
    font-size: .9vw;
    font-weight: bold;
  }

  .home__section2{
    display: flex;
    flex-direction: column-reverse;
    align-items: flex-end;
    margin-top: 2px;
  }

  .home__profile{
    background-color: #242F9B;
    width: 30px;
    height: 30px;
    border-radius: 50%;
    margin-right: 3vw;
  }
  .home__menu{
    display: flex;
    align-items: center;
    margin-right: 2vw;
    padding: 10px 15px 10px 15px;
  }
  .home__menu:hover{
    background-color: rgb(238, 238, 238);
  }
  .home__profile__text{
    /* margin-left: .9vw; */
    margin-top: 6px;
    text-align: center;
    color: white;
    font-weight: bold;
    font-size: 15px;
  }
  .profile__arrow{
    margin-left: -2.5vw;
    width:.5vw;
    height: 1vh;
  }
  .profile__setting{
    display: flex;
    background-color: white;
    flex-direction: column;
    padding: .5vw 1.5vw .5vw 1vw;
    margin-top: 10px;
    margin-right: 2vw;
    opacity: 0;
    z-index: 3;
  }
  .profile__setting.open {
    opacity: 1;
  }
  .profile__setting.open .profile__arrow{
    transform: rotate( 180deg );
  }
  .profile__func{
    margin-bottom: .7vh;
    display: flex;
    padding: .4vh 0 .4vh 0;
    align-items: center;
  }
  .profile__img{
    background-color: #242F9B;
    width: 1vw;
    height: 1vw;
    border-radius: 50%;

    margin-right: .5vw;
  }
  .profile__func__img{
    width: auto;
    height: 1.2vw;
    padding-right: .5vw;
  }
  .profile__func__txt{
    font-size:.7vw;
    font-weight: bold;
    cursor: default;
  }

  .sidebar.open ~ .app-contents{
    width: calc(100% - 13vw);
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
  .modal__profile__header{
    display: flex;
    justify-content: space-between;
    border-bottom: 1px white solid;
    padding-top: 1vw;
    padding-bottom: 15px;
    padding-left: 8px;
    margin: 15px 0px 0px 0px;
  }
  .modal__profile__header__account{
    margin-left: 10px;
  }
  .modal__profile__header__closeBtn{
    margin-right: 10px
  }
  .modal__profile__header__closeBtn:hover{
    cursor: default;
  }
  .modal__profile__main{
    display: flex;
  }
  .modal__profile__img{
    padding: 2.5vw;
    width: 10vw;
    height: auto;
    border-radius: 50%;
  }
  .modal__profile__setting{
    font-size: .8vw;
    padding-top: 5vh;
  }
  .setting__nickname,.setting__emailAddress, .setting__hashedPassword{
    margin: .6vh 0;
    margin-left: 1vw;
  }
  .setting__changeName, .setting__changeEmail, .setting__changePassword, .setting__delete{
    font-size: .1vw;
  }
  .setting__email, .setting__password{
    margin-top: 4vh;
  }
  .setting__delete{
    margin-top: 24vw;
    color: red;
  }
</style>

<template>
  <div>
    <section class = "sidebar" v-bind:class="{open: S_isActive}">
      <div class="sidebar__logo">
        <div class="sidebar__logo__name">
          <router-link to="/main">
            <img class="sidebar__logoImg" src="@/assets/logo_light.png" alt="">
          </router-link>
        </div>
        <i class="fa-solid fa-bars" @click="openSidebar"></i>
      </div>

      <ul class="sidebar__navList">
        <li class="i-sidebar__navList">
          <router-link to="/workflow">
            <div class="sidebar__mainMenu" @click="openMenu1">
            <!-- <i class="fa fa-home"></i> -->
              <img class ="mainMenu__icon" src="@/assets/workflow.png" alt="">
              <span class="sidebar__linksName">Workflows</span>
            </div>
          </router-link>
        </li>

        <li class="i-sidebar__navList">
          <router-link to="/files">
            <div class="sidebar__mainMenu" @click="openMenu2">
              <img class="mainMenu__icon" src="@/assets/files.png" alt="">
              <span class="sidebar__linksName">Files</span>
            </div>
          </router-link>
        </li>

        <li class="i-sidebar__navList">
          <div class="sidebar__mainMenu" @click="openMenu3">
            <!-- <i class="fa-solid fa-database mainMenu__icon"></i> -->
            <img class="mainMenu__icon" src="@/assets/database-file.png" alt="">
            <span class="sidebar__linksName">Data Sets</span>
          </div>
        </li>
      </ul>
    </section>

    <section class="home homeMain">
      <section class="home__section1">
        <template v-if="isUserLogin">
          <button @click="openSidebar">
            <i class="fa fa-solid fa-bars"></i>
          </button>
        </template>
        <button>
          <div class = "home__logo">
            <router-link to="/main">
              <img class="home__logoImg" src="@/assets/logo_dark.png" alt="">
            </router-link>
          </div>
        </button>

      </section>

      <section class="home__section2">
        <template v-if="isUserLogin">
          <ul class = "profile__setting" v-bind:class="{open: setting_isActive}">
            <li class="profile__func">
              <div class="profile__img">
                <!-- <div class="profile__text">J</div> -->
              </div>
              <div class="profile__func__txt"> {{ profile.email }} </div>
            </li>
            <li class="profile__func" @click="getProfile">
              <img class="profile__func__img" src="@/assets/person-profile-image.png" alt="">
              <div class="profile__func__txt"> internal profile </div>
            </li>
            <li class="profile__func">
              <img class="profile__func__img" src="@/assets/settings-gear.png" alt="">
              <div class="profile__func__txt"> settings </div>
            </li>
            <li class="profile__func">
              <img class="profile__func__img" src="@/assets/themes.png" alt="">
              <div class="profile__func__txt"> themes </div>
            </li>
            <li class="profile__func" style="border-top: 1px rgb(222,222,222) solid ">
              <img class="profile__func__img" src="@/assets/logout.png" alt="" style="padding-left:.2vw">
              <router-link to="/main">
                <div class="profile__func__txt" @click="logoutUser">
                  Logout
                </div>
              </router-link>
            </li>
          </ul>
          <div class="home__menu"  @click="openProfile">
            <div class="home__profile">
              <div class="home__profile__text"></div>
            </div>
            <img class = "profile__arrow" src="@/assets/arrow-bottom.png" alt="">
          </div>
        </template>

        <template v-else>
          <li class="home__login">
            <router-link to="/login" class="home__login__text">login</router-link>
          </li>
        </template>
        <Profile @click="modal=false" v-if="modal">
          <div slot="header" class="modal__profile__header">
            <div class="modal__profile__header__account">Account</div>
            <div class="modal__profile__header__closeBtn" @click="modal=false">X</div>
          </div>
          <div slot="body" class="modal__profile__main">
            <div>
              <img class="modal__profile__img" src="@/assets/no-profile-picture.png" alt="">
            </div>
            <div class="modal__profile__setting">
              <div class="setting__name">Name</div>
              <div class="setting__nickname">{{ profile.username }}</div>
              <div class="setting__changeName">Change name</div>
              <div class="setting__email">Email</div>
              <div class="setting__emailAddress">{{ profile.email }}</div>
              <div class="setting__changeEmail">Change mail</div>
              <div class="setting__password">Password</div>
              <div class="setting__hashedPassword">{{ profile.password }}</div>
              <div class="setting__changePassword">Change Password</div>
            </div>
            <div class="setting__delete">Delete account</div>
          </div>
        </Profile>
      </section>
    </section>
  </div>
</template>

<script>
import { getUser } from '@/api/index'
import { deleteCookie } from '@/utils/cookies'
import Profile from '@/components/profile'

export default {
  components: {Profile},
  data () {
    return {
      modal: false,
      S_isActive: false, // sidebar
      M1_isActive: false, // subMenu
      M2_isActive: false,
      M3_isActive: false,
      setting_isActive: false,
      profile: {
        username: null,
        email: null,
        password: '********'
      }
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
    async getProfile () {
      this.modal = true
      const userInfo = await getUser()
      this.$store.commit('setUserInfo', userInfo.data)
      console.log(userInfo.data)
      this.profile.username = userInfo.data.username
      this.profile.email = userInfo.data.email
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
    },
    openProfile () {
      this.setting_isActive = !this.setting_isActive
    }
  }
}

</script>

<style>
  @import '../../css/style.css';
</style>
