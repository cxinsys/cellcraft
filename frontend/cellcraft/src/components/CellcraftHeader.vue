<template>
   <section class="home-section" id="title">
      <section class="section1" style="margin-left: 1vw;">
        <button onclick="sidebar_open()">
          <i class="fa fa-solid fa-bars" id="btn"></i>
        </button>

        <button>
          <div id = "logo">LOGO</div>
        </button>
      </section>
      <section class="section2" style="display: flex;">

        <li class="menu">
          <a href="">
            <i class="fa fa-home"></i>
            <span class="links_name">Home</span>
          </a>
          <ul class="sub-menu">
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
        </li>
        <li class="menu">
          <a href="login.html">
            <span class="links_name">Login</span>
          </a>
        </li>
      </section>
    </section>
</template>

<script>
export default {

}
</script>

<style>
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
