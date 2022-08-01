<style scoped>
  @import '../../css/style.css';

  .profile {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 100vw;
    height: 100vh;
    position: relative;
    z-index: -1;
    background-color: rgba(0, 0, 0, 0.76);
  }
  .profile__page{
    width: 30vw;
    height: 30vw;
    background-color: rgb(0,29,110);
    position: relative;
    z-index: 100;
    color: white;
  }
  .profile__header{
    display: flex;
    justify-content: space-between;
    border-bottom: 1px white solid;
    /* margin-top: 2px; */
    padding-bottom: 15px;
    padding-left: 8px;
    margin: 15px 0px 0px 0px;
  }
  .profile__header__account{
    margin-left: 10px;
  }
  .profile__header__closeBtn{
    margin-right: 10px
  }
  .profile__main{
    display: flex;
  }
  .profile__img{
    padding: 2.5vw;
    width: 6vw;
    height: auto;
    border-radius: 50%;
  }
  .profile__setting{
    font-size: .8vw;
    padding-top: 5vh;
  }
  .setting__nickname,.setting__emailAddress{
    margin-top: .6vh;
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
  <section class="profile">
    <div class="profile__page">
      <div class="profile__header">
        <div class="profile__header__account">Account</div>
        <div class="profile__header__closeBtn">X</div>
      </div>
      <div class="profile__main">
        <div>
          <img class="profile__img" src="@/assets/no-profile-picture.png" alt="">
        </div>
        <div class="profile__setting">
          <div class="setting__name">Name</div>
          <div class="setting__nickname">Jisu Kang</div>
          <div class="setting__changeName">Change name</div>
          <div class="setting__email">Email</div>
          <div class="setting__emailAddress">email@uuu.com</div>
          <div class="setting__changeEmail">Change mail</div>
          <div class="setting__password">Password</div>
          <div class="setting__changePassword">Change Password</div>
        </div>
        <div class="setting__delete">Delete account</div>
      </div>
    </div>

 </section>
</template>

<script>
export default {

}
</script>
