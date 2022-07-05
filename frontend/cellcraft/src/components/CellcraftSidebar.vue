<template>
  <section class = "sidebar" id="sidebar">
    <div class="logo-details" style="display: flex;align-items: center; justify-content: space-between;">
      <div class="logo_name" style="margin-left: 1.5vw;">LOGO</div>
        <i class="fa-solid fa-bars" style="margin-right: .5vw;" id="btn"></i>
    </div>

    <ul class="nav-list">
      <li>
        <div class="mainMenu">
          <i class="fa fa-home"></i>
          <span class="links_name">Home</span>
        </div>
        <ul class="sub-menu">
          <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
          <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
          <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
        </ul>
      </li>

      <li>
        <div class="mainMenu">
          <i class="fa fa-user"></i>
          <span class="links_name">Content</span>
        </div>
        <ul class="sub-menu">
          <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
          <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
          <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
        </ul>
      </li>

      <li>
        <div class="mainMenu">
          <i class="fa fa-envelope-open icon"></i>
          <span class="links_name">Contact</span>
        </div>
        <ul class="sub-menu">
          <li><i class="fa-solid fa-heart"></i><a href="">Contact1</a></li>
          <li><i class="fa-solid fa-face-grin"></i><a href="">Contact2</a></li>
          <li><i class="fa-solid fa-heart"></i><a href="">Contact3</a></li>
        </ul>
      </li>

    </ul>
  </section>
</template>

<script>
export default {

}
</script>

<style>
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

</style>
