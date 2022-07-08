function saveAuthToCookie (value) {
  document.cookie = `scap_auth=${value}`
}

function getAuthFromCookie () {
  return document.cookie.replace(
    /(?:(?:^|.*;\s*)scap_auth\s*=\s*([^;]*).*$)|^.*$/,
    '$1'
  )
}

function deleteCookie () {
//   document.cookie = `${name}=; expires=Thu, 01 Jan 1970 00:00:01 GMT;`
  document.cookie = 'scap_auth=; expires=Thu, 01 Jan 1970 00:00:01 GMT;'
}

export {
  saveAuthToCookie,
  getAuthFromCookie,
  deleteCookie
}
