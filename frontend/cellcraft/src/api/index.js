import axios from 'axios'

axios.defaults.baseURL = 'http://127.0.0.1:8000' // 서버주소

// axios.defaults.headers.post['Content-Type'] = 'application/json;charset=utf-8'

// axios.defaults.headers.post['Access-Control-Allow-Origin'] = 'http://127.0.0.1:8000'

// axios.defaults.withCredentials = true

function registerUser (userData) {
  return axios.post('/routes/auth/register', userData)
}

function loginUser (userData) {
  return axios.post('/routes/auth/login/access-token', userData,
    {
      headers: {
        'Content-Type': 'multipart/form-data'
      }
    }
  )
}

export { registerUser, loginUser }
