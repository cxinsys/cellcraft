<template>
  <section class="container forms">
    <div class="form signup">
      <div class="form-content">
          <header>Sign Up</header>
          <form action="#">
            <div class="field input-field">
              <input type="email" placeholder="이메일을 입력해주세요" class="input">
            </div>

            <div class="field input-field">
              <input type="name" placeholder="이름을 입력해주세요" class="input">
            </div>

            <div class="field input-field">
              <input type="password" placeholder="비밀번호를 입력해주세요" class="password">
            </div>

            <div class="form-link">
              <p> 비밀번호 8자리 이상으로 설정해주세요</p>
            </div>

            <div class="field input-field">
              <input type="password" placeholder="비밀번호를 다시 입력해주세요" class="pwCheck">
            </div>

            <div class="field button-field">
              <button onclick="checkForm()">가입하기</button>
            </div>
          </form>

      </div>
    </div>
  </section>
</template>

<script>
export default {

}
</script>

<style>
  @import '../../css/style.css';
  *{
    box-sizing: border-box;
  }
  .container{
    height: 100vh;
    width: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    background-color: black;
    column-gap: 30px;
  }
  .form{
    position: absolute;
    max-width: 25vw;
    width: 100%;
    padding: 1.5vw;
    border-radius: .6vw;
    background: #FFF;
  }
  header{
    font-size: 1.5vw;
    font-weight: 600;
    text-align: center;
  }
  form{
    margin-top: 3vw;
  }
  .form .field{
    position: relative;
    height: 5vh;
    width: 100%;
    margin-top: 2vh;
    border-radius: .6vw;
  }
  .field input,
  .field button{
    height: 100%;
    width: 100%;
    border: none;
    font-size: .8vw;
    font-weight: 400;
    border-radius: .5vw;
  }
  .field input{
    outline: none;
    padding: 0 1vw;
    border: .1vw solid#CACACA;
  }
  .field input:focus{
    border-bottom-width: .2vw;
  }
  .field button{
    color: #fff;
    background-color: black;
    transition: all 0.3s ease;
    cursor: pointer;
  }
  /* .field button:hover{
    background-color: skyblue;
  } */
  .form-link{
    text-align: center;
    margin-top: 1vh;
  }
  .form-link a, .form-link p{
    font-size: .7vw;
    font-weight: 400;
    color: #232836;
  }

  .form a{
    color: #0171d3;
    text-decoration: none;
  }
  .form-content a:hover{
    text-decoration: underline;
  }
  .form-link span{
    font-size: .8vw;
    font-weight: 400;
    color: #232836;
  }
</style>
